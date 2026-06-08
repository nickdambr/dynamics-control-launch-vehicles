%% HM2 - Task 1: Powered Descent and Landing via Direct Collocation
%  Trapezoidal transcription, fixed-duration, minimum-fuel.
%  2D Cartesian, point-mass, no aerodynamics, flat Earth.
%
%  The OCP is solved in non-dimensional form (in the same spirit as HM1):
%  reference scales L_ref = y0, a_ref = g, t_ref = sqrt(L_ref/g),
%  V_ref = sqrt(g*L_ref), m_ref = m0, T_ref = m0*g.  The single residual
%  dimensionless parameter is V_c = V_ref/c (effective Tsiolkovsky number).
%  Internal solvers operate on the non-dim data; results are scaled back to
%  SI at the boundary for printing and plotting.
%
%  Reference: Homework 2 - Powered Descent Landing (Zavoli, April 2026)
%  Solver: fmincon (sqp).  No external dependency.

clear; close all; clc;

%% Problem data (Table 1, dimensional)
data.x0       = 1000;          % m
data.y0       = 3000;          % m
data.vx0      = 300;           % m/s
data.vy0      = -200;          % m/s
data.m0       = 2000;          % kg
data.g        = 9.81;          % m/s^2
data.Isp      = 225;           % s
data.g0       = 9.80665;       % m/s^2  (standard gravity)
data.c        = data.Isp * data.g0;
data.Tmin     = 0;             % N
data.Tmax     = 70000;         % N
data.theta_mx = deg2rad(60);   % glide-slope half-angle (radians, already nondim)

tf_nom = 38;                   % s
N      = 50;                   % collocation nodes

%% Non-dimensionalisation
[ref, dnd] = nondim(data);
fprintf('Non-dim reference scales:\n');
fprintf('  L_ref = %.1f m,  V_ref = %.2f m/s,  t_ref = %.3f s\n', ref.L, ref.V, ref.t);
fprintf('  m_ref = %.0f kg, T_ref = %.0f N,  V_c = V_ref/c = %.4f\n', ref.m, ref.T, dnd.Vc);

%% Sensitivity sweep on flight time (tf nominal +/- 5%)
tf_list = tf_nom * [0.95, 1.00, 1.05];
sols    = cell(numel(tf_list), 1);

for k = 1:numel(tf_list)
    fprintf('\n=== Solving for tf = %.2f s ===\n', tf_list(k));
    tf_nd  = tf_list(k) / ref.t;
    sol_nd = solve_trapcol(tf_nd, N, dnd);
    sols{k} = dim_sol(sol_nd, ref);
    fprintf('  final mass = %.2f kg, fuel used = %.2f kg\n', ...
        sols{k}.m_f, sols{k}.fuel);
end

%% Plots
plot_results(sols, tf_list, data);

%% Export figures
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

%% Summary table
fprintf('\n--- Sensitivity summary ---\n');
fprintf('%8s | %10s | %10s\n', 'tf [s]', 'm_f [kg]', 'fuel [kg]');
for k = 1:numel(tf_list)
    fprintf('%8.2f | %10.2f | %10.2f\n', tf_list(k), sols{k}.m_f, sols{k}.fuel);
end

%% =====================================================================
%  Helper functions
%  =====================================================================

function [ref, dnd] = nondim(d)
    % Reference scales (HM1-style choice: g and L set the units)
    ref.L = d.y0;                    % length    [m]
    ref.g = d.g;                     % accel.    [m/s^2]
    ref.t = sqrt(ref.L / ref.g);     % time      [s]
    ref.V = sqrt(ref.g * ref.L);     % velocity  [m/s]
    ref.m = d.m0;                    % mass      [kg]
    ref.T = ref.m * ref.g;           % thrust    [N]
    % Non-dim problem data
    dnd.x0       = d.x0  / ref.L;
    dnd.y0       = d.y0  / ref.L;
    dnd.vx0      = d.vx0 / ref.V;
    dnd.vy0      = d.vy0 / ref.V;
    dnd.m0       = d.m0  / ref.m;    % == 1
    dnd.Tmin     = d.Tmin / ref.T;
    dnd.Tmax     = d.Tmax / ref.T;
    dnd.Vc       = ref.V / d.c;      % the only residual nondim parameter
    dnd.theta_mx = d.theta_mx;
end

function sol = dim_sol(s_nd, ref)
    % Convert a non-dim solution struct back to SI units for output/plots.
    sol.t   = s_nd.t  * ref.t;
    sol.x   = s_nd.x  * ref.L;
    sol.y   = s_nd.y  * ref.L;
    sol.vx  = s_nd.vx * ref.V;
    sol.vy  = s_nd.vy * ref.V;
    sol.m   = s_nd.m  * ref.m;
    sol.Tx  = s_nd.Tx * ref.T;
    sol.Ty  = s_nd.Ty * ref.T;
    sol.Tmag = sqrt(sol.Tx.^2 + sol.Ty.^2);
    sol.tf   = s_nd.tf  * ref.t;
    sol.m_f  = s_nd.m_f * ref.m;
    sol.fuel = (s_nd.m0 - s_nd.m_f) * ref.m;
end

function sol = solve_trapcol(tf, N, d)
    % Trapezoidal direct collocation NLP (in non-dim variables).
    %   z = [x; y; vx; vy; m; Tx; Ty] stacked node-by-node, length 7*N.

    dt = tf / (N - 1);
    nz = 7 * N;
    idx = @(i) (i-1)*7 + (1:7);

    % --- Initial guess: linear state interpolation, hover thrust (T = m*1) ---
    z0 = zeros(nz, 1);
    for i = 1:N
        a  = (i-1) / (N-1);
        s  = idx(i);
        z0(s(1)) = (1-a) * d.x0;
        z0(s(2)) = (1-a) * d.y0;
        z0(s(3)) = (1-a) * d.vx0;
        z0(s(4)) = (1-a) * d.vy0;
        z0(s(5)) = d.m0 * (1 - 0.3*a);
        z0(s(6)) = 0;
        z0(s(7)) = d.m0;             % hover (gravity = 1 in nondim)
    end

    % --- Variable bounds ---
    lb = -inf(nz, 1);
    ub =  inf(nz, 1);
    for i = 1:N
        s = idx(i);
        lb(s(2)) = 0;
        lb(s(5)) = 1e-3;             % m strictly positive (small in nondim)
        ub(s(5)) = d.m0;
        lb(s(6)) = -d.Tmax;
        ub(s(6)) =  d.Tmax;
        lb(s(7)) = -d.Tmax;
        ub(s(7)) =  d.Tmax;
    end

    % --- Linear equality: boundary conditions ---
    s1 = idx(1);  sN = idx(N);
    Aeq = sparse(9, nz);
    beq = zeros(9, 1);
    Aeq(1, s1(1)) = 1;  beq(1) = d.x0;
    Aeq(2, s1(2)) = 1;  beq(2) = d.y0;
    Aeq(3, s1(3)) = 1;  beq(3) = d.vx0;
    Aeq(4, s1(4)) = 1;  beq(4) = d.vy0;
    Aeq(5, s1(5)) = 1;  beq(5) = d.m0;
    Aeq(6, sN(1)) = 1;  beq(6) = 0;
    Aeq(7, sN(2)) = 1;  beq(7) = 0;
    Aeq(8, sN(3)) = 1;  beq(8) = 0;
    Aeq(9, sN(4)) = 1;  beq(9) = 0;

    % --- Objective: maximize final mass ---
    iN_m  = (N-1)*7 + 5;
    f_obj = @(z) -z(iN_m);

    % --- Nonlinear constraints ---
    nlc = @(z) trap_nonlcon(z, N, dt, d);

    opts = optimoptions('fmincon', ...
        'Algorithm', 'sqp', ...
        'Display',   'final', ...
        'MaxIterations',          1000, ...
        'MaxFunctionEvaluations', 1e6, ...
        'OptimalityTolerance',    1e-5, ...
        'ConstraintTolerance',    1e-6, ...
        'StepTolerance',          1e-10);

    [z_opt, ~, exitflag] = fmincon(f_obj, z0, [], [], full(Aeq), beq, ...
                                   lb, ub, nlc, opts);

    if exitflag <= 0
        warning('fmincon did not converge cleanly (exitflag = %d)', exitflag);
    end

    % --- Unpack (still non-dim) ---
    Z = reshape(z_opt, 7, N).';
    sol.t   = linspace(0, tf, N).';
    sol.x   = Z(:,1);
    sol.y   = Z(:,2);
    sol.vx  = Z(:,3);
    sol.vy  = Z(:,4);
    sol.m   = Z(:,5);
    sol.Tx  = Z(:,6);
    sol.Ty  = Z(:,7);
    sol.tf  = tf;
    sol.m_f = sol.m(end);
    sol.m0  = d.m0;
end

function [c_ineq, c_eq] = trap_nonlcon(z, N, dt, d)
    Z = reshape(z, 7, N);

    % --- Dynamics RHS at every node (non-dim) ---
    f = zeros(5, N);
    for i = 1:N
        f(:,i) = dyn_rhs(Z(:,i), d.Vc);
    end

    % --- Trapezoidal defects ---
    defs = zeros(5, N-1);
    for k = 1:N-1
        defs(:,k) = Z(1:5, k+1) - Z(1:5, k) - 0.5*dt*(f(:,k) + f(:,k+1));
    end
    c_eq = defs(:);

    % --- Path constraints ---
    Tmag = sqrt(Z(6,:).^2 + Z(7,:).^2).';
    g_thr_lo = d.Tmin - Tmag;
    g_thr_hi = Tmag - d.Tmax;
    tt = tan(d.theta_mx);
    g_gs_pos = ( Z(1,:).' - tt*Z(2,:).');
    g_gs_neg = (-Z(1,:).' - tt*Z(2,:).');

    c_ineq = [g_thr_lo; g_thr_hi; g_gs_pos; g_gs_neg];
end

function dx = dyn_rhs(s, Vc)
    % Non-dim continuous dynamics.  State s = [x; y; vx; vy; m; Tx; Ty];
    % returns d/dt of [x; y; vx; vy; m].  Gravity is -1 and the only
    % remaining parameter is Vc = V_ref/c.
    Tmag = sqrt(s(6)^2 + s(7)^2);
    dx = [ s(3);
           s(4);
           s(6) / s(5);
           s(7) / s(5) - 1;
          -Vc * Tmag ];
end

function plot_results(sols, tf_list, d)
    colors = lines(numel(sols));
    lbl = arrayfun(@(t) sprintf('t_f = %.1f s', t), tf_list, 'UniformOutput', false);

    % --- Trajectory ---
    figure('Name','Trajectory','Position',[100 100 600 500]);
    hold on; grid on; axis equal;
    yy = linspace(0, max(cellfun(@(s) max(s.y), sols))*1.05, 50);
    xx = tan(d.theta_mx) * yy;
    plot( xx, yy, '--', 'Color',[0.5 0.5 0.5], 'HandleVisibility','off');
    plot(-xx, yy, '--', 'Color',[0.5 0.5 0.5], 'HandleVisibility','off');
    for k = 1:numel(sols)
        plot(sols{k}.x, sols{k}.y, '-', 'Color', colors(k,:), ...
            'LineWidth', 1.6, 'DisplayName', lbl{k});
    end
    plot(0, 0, 'k^', 'MarkerSize', 8, 'MarkerFaceColor','k', ...
        'HandleVisibility','off');
    xlabel('x  [m]'); ylabel('y  [m]');
    title('Powered descent trajectory  (dashed: glide-slope corridor)');
    legend('Location','best');

    % --- Thrust magnitude vs time ---
    figure('Name','Thrust magnitude','Position',[100 100 600 400]);
    hold on; grid on;
    for k = 1:numel(sols)
        plot(sols{k}.t, sols{k}.Tmag/1e3, '-', 'Color', colors(k,:), ...
            'LineWidth', 1.6, 'DisplayName', lbl{k});
    end
    yline(d.Tmax/1e3, 'k--', 'T_{max}', 'HandleVisibility','off');
    xlabel('t  [s]'); ylabel('|T|  [kN]');
    title('Thrust magnitude');
    legend('Location','best');

    % --- Mass vs time ---
    figure('Name','Mass','Position',[100 100 600 400]);
    hold on; grid on;
    for k = 1:numel(sols)
        plot(sols{k}.t, sols{k}.m, '-', 'Color', colors(k,:), ...
            'LineWidth', 1.6, 'DisplayName', lbl{k});
    end
    xlabel('t  [s]'); ylabel('m  [kg]');
    title('Vehicle mass');
    legend('Location','best');

    % --- Glide-slope angle vs time ---
    figure('Name','Glide-slope','Position',[100 100 600 400]);
    hold on; grid on;
    for k = 1:numel(sols)
        th = atan2(abs(sols{k}.x), sols{k}.y);
        plot(sols{k}.t, rad2deg(th), '-', 'Color', colors(k,:), ...
            'LineWidth', 1.6, 'DisplayName', lbl{k});
    end
    yline(rad2deg(d.theta_mx), 'k--', '\theta_{max}', 'HandleVisibility','off');
    xlabel('t  [s]'); ylabel('atan(|x|/y)  [deg]');
    title('Glide-slope angle');
    legend('Location','best');
end
