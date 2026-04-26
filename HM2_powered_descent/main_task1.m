%% HM2 - Task 1: Powered Descent and Landing via Direct Collocation
%  Trapezoidal transcription, fixed-duration, minimum-fuel.
%  2D Cartesian, point-mass, no aerodynamics, flat Earth.
%
%  Reference: Homework 2 - Powered Descent Landing (Zavoli, April 2026)
%  Solver: fmincon (sqp).  No external dependency.

clear; close all; clc;

%% Problem data (Table 1 of the assignment)
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
data.theta_mx = deg2rad(60);   % glide-slope half-angle

tf_nom = 38;                   % s
N      = 50;                   % collocation nodes

%% Sensitivity sweep on flight time (tf nominal +/- 5%)
tf_list = tf_nom * [0.95, 1.00, 1.05];
sols    = cell(numel(tf_list), 1);

for k = 1:numel(tf_list)
    fprintf('\n=== Solving for tf = %.2f s ===\n', tf_list(k));
    sols{k} = solve_trapcol(tf_list(k), N, data);
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

function sol = solve_trapcol(tf, N, d)
    % Trapezoidal direct collocation NLP.
    %   z = [x; y; vx; vy; m; Tx; Ty] stacked node-by-node, length 7*N.

    dt = tf / (N - 1);
    nz = 7 * N;

    % Index helper: idx(i) returns the 7-vector slice for node i
    idx = @(i) (i-1)*7 + (1:7);

    % --- Initial guess: linear interpolation, hover thrust ---
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
        z0(s(7)) = d.m0 * d.g;        % balance gravity at start
    end

    % --- Variable bounds ---
    lb = -inf(nz, 1);
    ub =  inf(nz, 1);
    for i = 1:N
        s = idx(i);
        lb(s(2)) = 0;                 % y >= 0  (cannot go underground)
        lb(s(5)) = 1;                 % m  > 0
        ub(s(5)) = d.m0;              % m <= m0
        lb(s(6)) = -d.Tmax;
        ub(s(6)) =  d.Tmax;
        lb(s(7)) = -d.Tmax;
        ub(s(7)) =  d.Tmax;
    end

    % --- Linear equality: boundary conditions ---
    s1 = idx(1);  sN = idx(N);
    Aeq = sparse(9, nz);
    beq = zeros(9, 1);
    Aeq(1, s1(1)) = 1;  beq(1) = d.x0;     % x(0)
    Aeq(2, s1(2)) = 1;  beq(2) = d.y0;     % y(0)
    Aeq(3, s1(3)) = 1;  beq(3) = d.vx0;    % vx(0)
    Aeq(4, s1(4)) = 1;  beq(4) = d.vy0;    % vy(0)
    Aeq(5, s1(5)) = 1;  beq(5) = d.m0;     % m(0)
    Aeq(6, sN(1)) = 1;  beq(6) = 0;        % x(tf)
    Aeq(7, sN(2)) = 1;  beq(7) = 0;        % y(tf)
    Aeq(8, sN(3)) = 1;  beq(8) = 0;        % vx(tf)
    Aeq(9, sN(4)) = 1;  beq(9) = 0;        % vy(tf)

    % --- Objective: minimize -m(tf)  (==  maximize final mass)  ---
    iN_m  = (N-1)*7 + 5;
    f_obj = @(z) -z(iN_m);

    % --- Nonlinear constraints: dynamics defects + thrust magnitude + glide-slope ---
    nlc = @(z) trap_nonlcon(z, N, dt, d, idx);

    opts = optimoptions('fmincon', ...
        'Algorithm', 'sqp', ...
        'Display',   'iter', ...
        'MaxIterations',          500, ...
        'MaxFunctionEvaluations', 2e5, ...
        'OptimalityTolerance',    1e-6, ...
        'ConstraintTolerance',    1e-6, ...
        'StepTolerance',          1e-10);

    [z_opt, ~, exitflag] = fmincon(f_obj, z0, [], [], Aeq, beq, lb, ub, nlc, opts);

    if exitflag <= 0
        warning('fmincon did not converge cleanly (exitflag = %d)', exitflag);
    end

    % --- Unpack ---
    Z = reshape(z_opt, 7, N).';     % N x 7
    sol.t   = linspace(0, tf, N).';
    sol.x   = Z(:,1);
    sol.y   = Z(:,2);
    sol.vx  = Z(:,3);
    sol.vy  = Z(:,4);
    sol.m   = Z(:,5);
    sol.Tx  = Z(:,6);
    sol.Ty  = Z(:,7);
    sol.Tmag= sqrt(sol.Tx.^2 + sol.Ty.^2);
    sol.tf  = tf;
    sol.m_f = sol.m(end);
    sol.fuel= d.m0 - sol.m_f;
end

function [c_ineq, c_eq] = trap_nonlcon(z, N, dt, d, idx)
    Z = reshape(z, 7, N);             % rows: x,y,vx,vy,m,Tx,Ty

    % --- Dynamics RHS at every node ---
    f = zeros(5, N);
    for i = 1:N
        f(:,i) = dyn_rhs(Z(:,i), d.c, d.g);
    end

    % --- Trapezoidal defects: x_{k+1} - x_k - dt/2 * (f_k + f_{k+1}) = 0 ---
    defs = zeros(5, N-1);
    for k = 1:N-1
        defs(:,k) = Z(1:5, k+1) - Z(1:5, k) - 0.5*dt*(f(:,k) + f(:,k+1));
    end
    c_eq = defs(:);

    % --- Thrust magnitude bounds: Tmin <= |T| <= Tmax ---
    Tmag = sqrt(Z(6,:).^2 + Z(7,:).^2).';
    g_thr_lo = d.Tmin - Tmag;         % <= 0
    g_thr_hi = Tmag - d.Tmax;         % <= 0

    % --- Glide-slope: |x| <= tan(theta_max) * y  ---
    tt = tan(d.theta_mx);
    g_gs_pos = ( Z(1,:).' - tt*Z(2,:).');   %  x - tt*y <= 0
    g_gs_neg = (-Z(1,:).' - tt*Z(2,:).');   % -x - tt*y <= 0

    c_ineq = [g_thr_lo; g_thr_hi; g_gs_pos; g_gs_neg];
end

function dx = dyn_rhs(s, c, g)
    % State s = [x; y; vx; vy; m; Tx; Ty]; returns d/dt of [x; y; vx; vy; m].
    Tmag = sqrt(s(6)^2 + s(7)^2);
    dx = [ s(3);
           s(4);
           s(6) / s(5);
           s(7) / s(5) - g;
          -Tmag / c ];
end

function plot_results(sols, tf_list, d)
    colors = lines(numel(sols));
    lbl = arrayfun(@(t) sprintf('t_f = %.1f s', t), tf_list, 'UniformOutput', false);

    % --- Trajectory ---
    figure('Name','Trajectory','Position',[100 100 600 500]);
    hold on; grid on; axis equal;
    % glide-slope corridor (drawn from origin)
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
