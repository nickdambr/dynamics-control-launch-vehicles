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

% Initial guess (improved formulation: 4 unknowns, costates normalized to lam_m0=1)
z_guess = [0.6; 3.8; 14; 0.30];

fprintf('\nPhase 2: Solving BVP for optimal burn...\n');
[z_sol, ~, ef] = fsolve(@(z) shooting2(z, p, opts_ode), z_guess, opts_fs);

if ef <= 0
    fprintf('WARNING: fsolve did not converge. Trying alternative guess...\n');
    z_guess2 = [0.4; 3.0; 10; 0.35];
    [z_sol, ~, ef] = fsolve(@(z) shooting2(z, p, opts_ode), z_guess2, opts_fs);
end

if ef > 0
    % Extract solution
    pp = p;
    pp.lam_vx0 = z_sol(1);
    pp.lam_vy0 = z_sol(2);
    pp.lam_y   = z_sol(3);
    t_burn = z_sol(4);

    ic2 = [p.x0; p.y0; p.vx0; p.vy0; p.m0; 1];
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
    ax_main = axes;
    plot(ax_main, x_vert, y_vert, 'r-', 'LineWidth', 2, 'DisplayName', 'Vertical climb');
    hold(ax_main, 'on');
    plot(ax_main, x_burn, y_burn, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Optimal burn');
    xlabel(ax_main, 'x (nondim)'); ylabel(ax_main, 'y (nondim)');
    title(ax_main, sprintf('Task 2: Trajectory (Q = %.1f, y_f = %.2f)', Q, yf));
    legend(ax_main, 'Location', 'best'); grid(ax_main, 'on');

    % --- Inset: zoom on the vertical-climb + pitch-over knee ---
    % The vertical climb (x = 0, 0 <= y <= y1) spans y1 = 1e-4, i.e. 1/400 of the
    % full ascent, so it is invisible at full scale; an inset magnifies it.
    y_zoom = 5 * y1;                                    % tight zoom on the climb
    kk_in  = find(y_burn <= y_zoom, 1, 'last');
    if isempty(kk_in) || kk_in < 2, kk_in = min(40, numel(x_burn)); end
    x_zmax = max(x_burn(1:kk_in)) * 1.10 + 1e-6;

    % dashed rectangle on the main axes marking the magnified region
    rectangle(ax_main, 'Position', [-0.02*x_zmax, 0, 1.04*x_zmax, y_zoom], ...
              'EdgeColor', [0.4 0.4 0.4], 'LineStyle', '--', 'LineWidth', 0.8);

    ax_in = axes('Position', [0.15 0.55 0.27 0.32]);    % inset (norm. figure units)
    plot(ax_in, x_vert, y_vert, 'r-', 'LineWidth', 2); hold(ax_in, 'on');
    plot(ax_in, x_burn, y_burn, 'b-', 'LineWidth', 1.5);
    plot(ax_in, 0, y_1, 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k');  % climb->burn handover
    xlim(ax_in, [-0.02*x_zmax, x_zmax]); ylim(ax_in, [0, y_zoom]);
    grid(ax_in, 'on'); box(ax_in, 'on'); set(ax_in, 'FontSize', 8);
    title(ax_in, 'vertical-climb zoom', 'FontSize', 8);

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
    % Vertical phase: thrust phi = 90 deg (solid, so the blue curve starts at the
    % climb); velocity psi = 90 deg too (the vehicle flies straight up)
    plot([0 t1], [90 90], 'b-',  'LineWidth', 1.5, 'HandleVisibility','off');
    plot([0 t1], [90 90], 'r--', 'LineWidth', 1.5, 'HandleVisibility','off');
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
% Shooting residual for the burn arc after the vertical climb, in the improved
% formulation (lam_m0 = 1 normalization; H = 0 imposed algebraically at the
% START of the burn, where the state (0, y1, 0, vy1, m1) is known).
%   z0 = [lam_vx0; lam_vy0; lam_y; t_burn]

    lam_vx0 = z0(1);
    lam_vy0 = z0(2);
    lam_y   = z0(3);
    t_burn  = z0(4);

    if t_burn <= 0 || t_burn > 2
        res = 1e6 * ones(4,1);
        return;
    end

    pp = p;
    pp.lam_vx0 = lam_vx0;
    pp.lam_vy0 = lam_vy0;
    pp.lam_y   = lam_y;

    ic = [p.x0; p.y0; p.vx0; p.vy0; p.m0; 1];   % lam_m0 = 1 (normalization)

    try
        [~, Z] = ode45(@(t,z) ode_burn(t, z, pp), [0 t_burn], ic, opts_ode);
        zf = Z(end,:);
    catch
        res = 1e6 * ones(4,1);
        return;
    end

    % H = 0 at the start of the burn (state known): vx0 = 0, vy0 = vy1,
    % m0 = m1, lam_m0 = 1, lam_x = 0.
    Lam0 = sqrt(lam_vx0^2 + lam_vy0^2);
    H0 = lam_y * p.vy0 + (p.T / p.m0) * Lam0 - lam_vy0 - p.T / p.c;

    res = [zf(2) - p.yf;
           zf(3) - 1;
           zf(4);
           H0];
end
