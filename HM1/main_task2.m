%% TASK 2 - Vertical climb + optimal burn arc
%  Vertical climb to y1 = 0.0001, then optimal gravity turn
%  Compare payload with Task 1 (no vertical climb)

clear; close all; clc;

%% Parameters
c   = 0.6;
eta = 0.1;
y1  = 0.0001;   % vertical climb altitude
yf  = 0.04;     % target orbit altitude
Q   = 2;         % reference mass flow rate (can be changed)

opts_ode = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
opts_fs  = optimoptions('fsolve', 'Display', 'iter', ...
    'MaxIterations', 500, 'MaxFunctionEvaluations', 5000, ...
    'FunctionTolerance', 1e-10, 'StepTolerance', 1e-10);

%% Phase 1: Vertical climb (phi = pi/2) until y = y1
T = c * Q;
fprintf('Phase 1: Vertical climb to y1 = %.4f\n', y1);
fprintf('  T/W = %.3f\n', T);

% Integrate vertical ODE: z = [y; vy; m]
ic_vert = [0; 0; 1];
opts_vert = odeset('RelTol', 1e-12, 'AbsTol', 1e-14, ...
    'Events', @(t,z) event_altitude(t, z, y1));
[T_vert, Z_vert] = ode45(@(t,z) ode_vertical(t, z, T, Q), [0 1], ic_vert, opts_vert);

t1  = T_vert(end);
y_1 = Z_vert(end,1);
vy1 = Z_vert(end,2);
m1  = Z_vert(end,3);

fprintf('  t1 = %.6f\n', t1);
fprintf('  vy(t1) = %.6f\n', vy1);
fprintf('  m(t1)  = %.6f\n', m1);

%% Phase 2: Optimal burn (same BVP as Task 1 but from vertical climb state)
% State at start of Phase 2: x=0, y=y1, vx=0, vy=vy1, m=m1
% BVP unknowns: [lam_vx0, lam_vy0, lam_y, lam_m0, t_burn]
% BVP targets:  [y(tf)-yf, vx(tf)-1, vy(tf), lam_m(tf)-1, H(tf)]

p.c  = c;
p.Q  = Q;
p.T  = T;
p.yf = yf;
p.x0  = 0;
p.y0  = y_1;
p.vx0 = 0;
p.vy0 = vy1;
p.m0  = m1;

% Initial guess (similar to Task 1, adjusted for vertical climb IC)
z_guess = [0.22; 1.3; 4.5; -0.5; 0.30];

fprintf('\nPhase 2: Solving BVP for optimal burn...\n');
[z_sol, ~, ef] = fsolve(@(z) shooting2(z, p, opts_ode), z_guess, opts_fs);

if ef <= 0
    fprintf('WARNING: fsolve did not converge. Trying alternative guess...\n');
    z_guess2 = [0.15; 1.0; 3.5; -0.3; 0.35];
    [z_sol, ~, ef] = fsolve(@(z) shooting2(z, p, opts_ode), z_guess2, opts_fs);
end

if ef > 0
    % Extract solution
    pp = p;
    pp.lam_vx0 = z_sol(1);
    pp.lam_vy0 = z_sol(2);
    pp.lam_y   = z_sol(3);
    t_burn = z_sol(5);

    ic2 = [p.x0; p.y0; p.vx0; p.vy0; p.m0; z_sol(4)];
    [T2, Z2] = ode45(@(t,z) ode_burn(t, z, pp), linspace(0, t_burn, 500), ic2, opts_ode);

    mf_task2 = Z2(end, 5);
    tf_total = t1 + t_burn;

    fprintf('\n=== Task 2 Results (Q = %.1f, yf = %.2f) ===\n', Q, yf);
    fprintf('  Vertical climb time: %.6f\n', t1);
    fprintf('  Burn time:           %.6f\n', t_burn);
    fprintf('  Total flight time:   %.6f\n', tf_total);
    fprintf('  Final mass:          %.6f\n', mf_task2);

    % Payload comparison
    payload_task2 = mf_task2 - eta * (1 - mf_task2) / (1 + eta);
    fprintf('  Payload (approx):    %.6f\n', mf_task2 * (1 + eta) - eta);

    %% Plots
    % Combine vertical and burn trajectories
    x_vert = zeros(size(T_vert));
    y_vert = Z_vert(:,1);

    x_burn = Z2(:,1);
    y_burn = Z2(:,2);
    vx_burn = Z2(:,3);
    vy_burn = Z2(:,4);

    figure('Name','Task 2 - Trajectory');
    plot(x_vert, y_vert, 'r-', 'LineWidth', 2, 'DisplayName', 'Vertical climb');
    hold on;
    plot(x_burn, y_burn, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Optimal burn');
    xlabel('x (nondim)'); ylabel('y (nondim)');
    title(sprintf('Task 2: Trajectory (Q = %.1f, y_f = %.2f)', Q, yf));
    legend('Location','best'); grid on;

    % Angles
    phi_burn = zeros(size(T2));
    psi_burn = zeros(size(T2));
    for kk = 1:length(T2)
        lam_vy_k = z_sol(2) - z_sol(3) * T2(kk);
        phi_burn(kk) = atan2(lam_vy_k, z_sol(1));
        V_k = sqrt(vx_burn(kk)^2 + vy_burn(kk)^2);
        if V_k > 1e-10
            psi_burn(kk) = atan2(vy_burn(kk), vx_burn(kk));
        else
            psi_burn(kk) = phi_burn(kk);
        end
    end

    figure('Name','Task 2 - Angles');
    hold on; grid on;
    % Vertical phase: phi = 90 deg
    plot([0 t1], [90 90], 'b--', 'LineWidth', 1, 'HandleVisibility','off');
    plot([0 t1], [90 90], 'r--', 'LineWidth', 1, 'HandleVisibility','off');
    % Burn phase
    plot(t1 + T2, rad2deg(phi_burn), 'b-', 'LineWidth', 1.5, 'DisplayName', '\phi (thrust)');
    plot(t1 + T2, rad2deg(psi_burn), 'r--', 'LineWidth', 1.5, 'DisplayName', '\psi (velocity)');
    xlabel('Time (nondim)'); ylabel('Angle (deg)');
    title('Task 2: Thrust and velocity angles');
    legend('Location','best');
else
    fprintf('ERROR: BVP did not converge.\n');
end

%% ===================== EXPORT FIGURES =====================
fig_dir = fullfile(fileparts(mfilename('fullpath')), 'figures');
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end
slugify = @(s) lower(regexprep(s, '[^a-zA-Z0-9]+', '_'));
fig_handles = findobj(groot, 'Type', 'figure');
for kk = 1:numel(fig_handles)
    nm = get(fig_handles(kk), 'Name');
    if isempty(nm); nm = sprintf('fig%d', kk); end
    try
        theme(fig_handles(kk), 'light');    % force light theme (ignore desktop dark mode)
        drawnow;
    catch
        set(fig_handles(kk), 'Color', 'w'); % fallback for pre-R2025a MATLAB
    end
    exportgraphics(fig_handles(kk), ...
        fullfile(fig_dir, ['task2_' slugify(nm) '.png']), 'Resolution', 200);
end

%% ===================== LOCAL FUNCTIONS =====================

function dz = ode_vertical(t, z, T, Q)
% Vertical climb ODE: z = [y; vy; m], phi = pi/2
    vy = z(2); m = z(3);
    dz = [vy; T/m - 1; -Q];
end

function [value, isterminal, direction] = event_altitude(t, z, y_target)
    value = z(1) - y_target;
    isterminal = 1;
    direction = 1;
end

function res = shooting2(z0, p, opts_ode)
% Shooting residual for burn arc after vertical climb
%   z0 = [lam_vx0; lam_vy0; lam_y; lam_m0; t_burn]

    lam_vx0 = z0(1);
    lam_vy0 = z0(2);
    lam_y   = z0(3);
    lam_m0  = z0(4);
    t_burn  = z0(5);

    if t_burn <= 0 || t_burn > 2
        res = 1e6 * ones(5,1);
        return;
    end

    pp = p;
    pp.lam_vx0 = lam_vx0;
    pp.lam_vy0 = lam_vy0;
    pp.lam_y   = lam_y;

    ic = [p.x0; p.y0; p.vx0; p.vy0; p.m0; lam_m0];

    try
        [~, Z] = ode45(@(t,z) ode_burn(t, z, pp), [0 t_burn], ic, opts_ode);
        zf = Z(end,:);
    catch
        res = 1e6 * ones(5,1);
        return;
    end

    y_f    = zf(2);
    vx_f   = zf(3);
    vy_f   = zf(4);
    m_f    = zf(5);
    lam_mf = zf(6);

    lam_vy_f   = lam_vy0 - lam_y * t_burn;
    lam_v_norm = sqrt(lam_vx0^2 + lam_vy_f^2);

    H_f = lam_y * vy_f + (p.T / m_f) * lam_v_norm - lam_vy_f - p.Q * lam_mf;

    res = [y_f - p.yf;
           vx_f - 1;
           vy_f;
           lam_mf - 1;
           H_f];
end
