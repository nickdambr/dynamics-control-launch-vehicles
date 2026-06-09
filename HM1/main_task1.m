%% TASK 1 - Single burn arc: minimum-fuel ascent trajectory
%  (a) mf vs Q for yf = 0.04, 0.05, 0.06
%  (b) Velocity losses vs Q for yf = 0.04
%  (c) Trajectory and angles for optimal Q at yf = 0.04

clear; close all; clc;

%% Parameters
c   = 0.6;          % effective exhaust velocity (nondim)
eta = 0.1;          % structural coefficient ms/mp

yf_vec = [0.04, 0.05, 0.06];
Q_vec  = linspace(1.8, 7, 80);   % Q > 1/c ~ 1.667 for liftoff

opts_ode = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
opts_fs  = optimoptions('fsolve', 'Display', 'off', ...
    'MaxIterations', 500, 'MaxFunctionEvaluations', 5000, ...
    'FunctionTolerance', 1e-10, 'StepTolerance', 1e-10);

%% Storage
nQ  = length(Q_vec);
nyf = length(yf_vec);
mf_results  = nan(nQ, nyf);
sol_results = cell(nQ, nyf);

%% Solve BVP for each yf, sweeping Q with continuation
for jj = 1:nyf
    yf = yf_vec(jj);
    fprintf('--- Solving for yf = %.2f ---\n', yf);

    % Setup parameter struct
    p.c  = c;
    p.yf = yf;

    % Start from Q near 3 (good T/W ratio)
    [~, idx0] = min(abs(Q_vec - 3));

    % Initial guess: [lam_vx0, lam_vy0, lam_y, lam_m0, tf]
    if jj == 1
        z_guess = [0.6; 3.8; 14; 0.30];
    else
        % Use solution from previous yf at same Q
        z_guess = sol_results{idx0, jj-1};
        if isempty(z_guess)
            z_guess = [0.6; 3.8; 14; 0.30];
        end
    end

    % Solve at starting Q
    p.Q = Q_vec(idx0);
    p.T = c * p.Q;
    [z_sol, ~, ef] = fsolve(@(z) shooting1(z, p, opts_ode), z_guess, opts_fs);

    if ef > 0
        [~, Z] = ode45(@(t,z) ode_burn(t,z,set_costates(p,z_sol)), ...
                        [0 z_sol(4)], [0;0;0;0;1;1], opts_ode);
        mf_results(idx0,jj) = Z(end,5);
        sol_results{idx0,jj} = z_sol;
        fprintf('  Q=%.2f  mf=%.5f  converged\n', Q_vec(idx0), Z(end,5));
    else
        fprintf('  Q=%.2f  FAILED at starting point\n', Q_vec(idx0));
        continue;
    end

    % Sweep forward (increasing Q)
    z_prev = z_sol;
    for ii = idx0+1:nQ
        p.Q = Q_vec(ii);
        p.T = c * p.Q;
        [z_sol, ~, ef] = fsolve(@(z) shooting1(z, p, opts_ode), z_prev, opts_fs);
        if ef > 0
            [~, Z] = ode45(@(t,z) ode_burn(t,z,set_costates(p,z_sol)), ...
                            [0 z_sol(4)], [0;0;0;0;1;1], opts_ode);
            mf_results(ii,jj) = Z(end,5);
            sol_results{ii,jj} = z_sol;
            z_prev = z_sol;
        end
    end

    % Sweep backward (decreasing Q)
    z_prev = sol_results{idx0,jj};
    for ii = idx0-1:-1:1
        p.Q = Q_vec(ii);
        p.T = c * p.Q;
        [z_sol, ~, ef] = fsolve(@(z) shooting1(z, p, opts_ode), z_prev, opts_fs);
        if ef > 0
            [~, Z] = ode45(@(t,z) ode_burn(t,z,set_costates(p,z_sol)), ...
                            [0 z_sol(4)], [0;0;0;0;1;1], opts_ode);
            mf_results(ii,jj) = Z(end,5);
            sol_results{ii,jj} = z_sol;
            z_prev = z_sol;
        end
    end
end

%% ===== Optimal Q* per target altitude =====
fprintf('\n=== Optimal mass-flow rate Q* per target altitude ===\n');
Q_star  = nan(1, nyf);
mf_star = nan(1, nyf);
for jj = 1:nyf
    [mf_star(jj), idx] = max(mf_results(:, jj));
    Q_star(jj) = Q_vec(idx);
    fprintf('  yf=%.2f:  Q*=%.4f  mf*=%.5f  payload*=%.5f\n', ...
        yf_vec(jj), Q_star(jj), mf_star(jj), mf_star(jj)*(1+eta)-eta);
end

%% ===== PLOT 1a: Final mass vs Q =====
figure('Name','Task 1a - Final mass vs Q');
hold on; grid on;
colors = {'b','r','k'};
for jj = 1:nyf
    valid = ~isnan(mf_results(:,jj));
    plot(Q_vec(valid), mf_results(valid,jj), [colors{jj} '-o'], ...
        'MarkerSize', 3, 'DisplayName', sprintf('y_f = %.2f', yf_vec(jj)));
end
xlabel('Mass flow rate Q'); ylabel('Final mass m_f');
title('Task 1a: Final mass vs mass flow rate');
legend('Location','best');

%% ===== PLOT 1b: Velocity losses vs Q for yf = 0.04 =====
jj_ref = 1;  % yf = 0.04
yf = yf_vec(jj_ref);
Wd_vec = nan(nQ,1);
Wg_vec = nan(nQ,1);
Wt_vec = nan(nQ,1);

for ii = 1:nQ
    z_sol = sol_results{ii, jj_ref};
    if isempty(z_sol), continue; end

    p.c = c; p.yf = yf;
    p.Q = Q_vec(ii); p.T = c * p.Q;
    pp = set_costates(p, z_sol);
    tf = z_sol(4);

    % Integrate with losses
    ic8 = [0; 0; 0; 0; 1; 1; 0; 0];
    [~, Z8] = ode45(@(t,z) ode_burn_losses(t, z, pp), [0 tf], ic8, opts_ode);
    Wd_vec(ii) = Z8(end,7);
    Wg_vec(ii) = Z8(end,8);
    mf_ii = Z8(end,5);
    Wt_vec(ii) = c * log(1/mf_ii) - 1;  % total = DV_ideal - V_final
end

figure('Name','Task 1b - Velocity losses');
hold on; grid on;
valid = ~isnan(Wd_vec);
plot(Q_vec(valid), Wd_vec(valid), 'b-o', 'MarkerSize', 3, 'DisplayName', 'W_d (misalignment)');
plot(Q_vec(valid), Wg_vec(valid), 'r-o', 'MarkerSize', 3, 'DisplayName', 'W_g (gravity)');
plot(Q_vec(valid), Wt_vec(valid), 'k-o', 'MarkerSize', 3, 'DisplayName', 'W_{tot}');
xlabel('Mass flow rate Q'); ylabel('Velocity loss (nondim)');
title(sprintf('Task 1b: Velocity losses (y_f = %.2f)', yf));
legend('Location','best');

%% ===== PLOT 1c: Trajectory and angles for optimal Q at yf = 0.04 =====
% Find Q that maximizes mf
[mf_opt, idx_opt] = max(mf_results(:, jj_ref));
Q_opt = Q_vec(idx_opt);
z_opt = sol_results{idx_opt, jj_ref};
fprintf('\nOptimal Q = %.3f  ->  mf = %.5f  (yf = %.2f)\n', Q_opt, mf_opt, yf);

p.c = c; p.yf = yf;
p.Q = Q_opt; p.T = c * Q_opt;
pp = set_costates(p, z_opt);
tf = z_opt(4);

% Dense integration for plotting
ic8 = [0; 0; 0; 0; 1; 1; 0; 0];
[T_sol, Z_sol] = ode45(@(t,z) ode_burn_losses(t, z, pp), linspace(0, tf, 500), ic8, opts_ode);

x_traj  = Z_sol(:,1);
y_traj  = Z_sol(:,2);
vx_traj = Z_sol(:,3);
vy_traj = Z_sol(:,4);
m_traj  = Z_sol(:,5);

% Compute angles
phi_traj = zeros(size(T_sol));
psi_traj = zeros(size(T_sol));
for kk = 1:length(T_sol)
    lam_vy_k = z_opt(2) - z_opt(3) * T_sol(kk);
    phi_traj(kk) = atan2(lam_vy_k, z_opt(1));
    V_k = sqrt(vx_traj(kk)^2 + vy_traj(kk)^2);
    if V_k > 1e-10
        psi_traj(kk) = atan2(vy_traj(kk), vx_traj(kk));
    else
        psi_traj(kk) = phi_traj(kk);
    end
end

figure('Name','Task 1c - Trajectory');
plot(x_traj, y_traj, 'b-', 'LineWidth', 1.5);
xlabel('x (nondim)'); ylabel('y (nondim)');
title(sprintf('Task 1c: Trajectory (Q = %.2f, y_f = %.2f)', Q_opt, yf));
grid on; axis equal;

figure('Name','Task 1c - Thrust and velocity angles');
hold on; grid on;
plot(T_sol, rad2deg(phi_traj), 'b-', 'LineWidth', 1.5, 'DisplayName', '\phi (thrust)');
plot(T_sol, rad2deg(psi_traj), 'r--', 'LineWidth', 1.5, 'DisplayName', '\psi (velocity)');
xlabel('Time (nondim)'); ylabel('Angle (deg)');
title(sprintf('Task 1c: Thrust and velocity angles (Q = %.2f)', Q_opt));
legend('Location','best');

fprintf('\nTask 1 complete.\n');

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
        fullfile(fig_dir, ['task1_' slugify(nm) '.png']), 'Resolution', 200);
end

%% ===================== LOCAL FUNCTIONS =====================

function res = shooting1(z0, p, opts_ode)
% Shooting residual for a single burn arc, in the improved formulation of the
% course notes (HW1, "migliorie numeriche"):
%   (i)  costates are NORMALIZED by fixing lam_m0 = 1 (H is homogeneous of
%        degree 1 in lambda, so the solution is defined up to a scale);
%   (ii) the free-time condition H = 0 is imposed at the INITIAL instant,
%        where it is purely algebraic (vx0 = vy0 = 0) and accumulates no
%        integration error, instead of at t_f.
% Each move removes one unknown/condition, leaving 4 unknowns and 4 residuals;
% lam_m never enters the residual.
%   z0 = [lam_vx0; lam_vy0; lam_y; tf]

    lam_vx0 = z0(1);
    lam_vy0 = z0(2);
    lam_y   = z0(3);
    tf      = z0(4);

    if tf <= 0 || tf > 2
        res = 1e6 * ones(4,1);
        return;
    end

    pp = p;
    pp.lam_vx0 = lam_vx0;
    pp.lam_vy0 = lam_vy0;
    pp.lam_y   = lam_y;

    ic = [0; 0; 0; 0; 1; 1];   % state + lam_m0 = 1 (normalization)

    try
        [~, Z] = ode45(@(t,z) ode_burn(t, z, pp), [0 tf], ic, opts_ode);
        zf = Z(end,:);
    catch
        res = 1e6 * ones(4,1);
        return;
    end

    % H = 0 at t0 (autonomous => H const = 0). With vx0 = vy0 = 0, m0 = 1 and
    % lam_m0 = 1:   H0 = -lam_vy0 + T*( |lam_v0| - 1/c ).
    Lam0 = sqrt(lam_vx0^2 + lam_vy0^2);
    H0   = -lam_vy0 + p.T * (Lam0 - 1/p.c);

    res = [zf(2) - p.yf;    % y(tf)  = yf
           zf(3) - 1;       % vx(tf) = 1
           zf(4);           % vy(tf) = 0
           H0];             % H(0)   = 0  (free final time)
end

function pp = set_costates(p, z_sol)
% Pack costate parameters into struct
    pp = p;
    pp.lam_vx0 = z_sol(1);
    pp.lam_vy0 = z_sol(2);
    pp.lam_y   = z_sol(3);
end

function dz = ode_burn_losses(t, z, p)
% Extended ODE with velocity losses
%   z = [x; y; vx; vy; m; lam_m; Wd; Wg]

    vx = z(3); vy = z(4); m = z(5);

    lam_vx = p.lam_vx0;
    lam_vy = p.lam_vy0 - p.lam_y * t;
    lam_v_norm = sqrt(lam_vx^2 + lam_vy^2);

    phi = atan2(lam_vy, lam_vx);

    V = sqrt(vx^2 + vy^2);
    if V > 1e-12
        psi = atan2(vy, vx);
    else
        psi = phi;
    end

    dz = zeros(8,1);
    dz(1) = vx;
    dz(2) = vy;
    dz(3) = (p.T / m) * cos(phi);
    dz(4) = (p.T / m) * sin(phi) - 1;
    dz(5) = -p.Q;
    dz(6) = (p.T / m^2) * lam_v_norm;
    dz(7) = (p.T / m) * (1 - cos(phi - psi));   % dWd/dt
    dz(8) = sin(psi);                             % dWg/dt
end
