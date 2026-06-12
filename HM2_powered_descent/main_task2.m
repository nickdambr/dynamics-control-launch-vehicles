%% HM2 - Task 2: ZOH discretization (Optional)
%  Two parallel implementations of the ZOH-based transcription:
%    (a) Nonlinear ZOH with RK4 propagation
%        -- multiple-shooting NLP with x_{k+1} = RK4(x_k, u_k, dt)
%    (b) LTV-linearized ZOH per Appendix A, wrapped in an SCvx outer loop
%        -- linearize about a reference, compute discrete-time matrices via
%           the auxiliary ODE of Appendix A, solve the LTV NLP, update the
%           reference, repeat to convergence.
%
%  Both share the trapezoidal Task 1 baseline as a reference point and are
%  solved in **non-dimensional form** (HM1-style scaling: L_ref = y0,
%  a_ref = g, t_ref = sqrt(L_ref/g), V_ref = sqrt(g*L_ref), m_ref = m0,
%  T_ref = m0*g, with V_c = V_ref/c the only residual nondim parameter).
%  Solutions are scaled back to SI for printing and plotting.
%
%  Validation: forward-integrate the optimized ZOH controls with ode45 and
%  compare against each transcription's discretized trajectory.
%
%  Reference: Homework 2 - Powered Descent Landing (Zavoli, April 2026)
%  Solver: fmincon (sqp).  No external dependency.

clear; close all; clc;

%% Problem data (Table 1, dimensional)
data.x0       = 1000;     data.y0       = 3000;
data.vx0      = 300;      data.vy0      = -200;
data.m0       = 2000;
data.g        = 9.81;
data.Isp      = 225;      data.g0       = 9.80665;
data.c        = data.Isp * data.g0;
data.Tmin     = 0;        data.Tmax     = 70000;
data.theta_mx = deg2rad(60);

tf       = 38;     % flight time (fixed) [s]
N        = 50;     % grid nodes
n_sub    = 2;      % RK4 substeps per ZOH interval (variant a)
scvx_max = 15;     % SCvx maximum iterations    (variant b)
scvx_tol = 1e-3;   % SCvx convergence tolerance (variant b, nondim)

%% Non-dimensionalisation
[ref, dnd] = nondim(data);
fprintf('Non-dim reference scales:\n');
fprintf('  L_ref = %.1f m,  V_ref = %.2f m/s,  t_ref = %.3f s\n', ref.L, ref.V, ref.t);
fprintf('  m_ref = %.0f kg, T_ref = %.0f N,  V_c = V_ref/c = %.4f\n', ref.m, ref.T, dnd.Vc);
tf_nd = tf / ref.t;

% SCvx trust-region radii (in non-dim units)
trust = struct('pos', 0.17, 'vel', 0.6, 'mass', 0.1, 'thrust', 1.0);

%% YALMIP + ECOS availability check
yalmip_ok = exist('yalmip', 'file') && exist('ecos', 'file');
if ~yalmip_ok
    warning('YALMIP or ECOS not found — YALMIP/ECOS variant will be skipped.');
end

%% Solve all three transcriptions
fprintf('\n=== Trapezoidal transcription (Task 1 baseline) ===\n');
t0 = tic;
sol_trap_nd = solve_trap(tf_nd, N, dnd);
t_trap = toc(t0);
sol_trap = dim_sol(sol_trap_nd, ref);
fprintf('  m_f = %.2f kg, fuel = %.2f kg, wall = %.1f s\n', ...
    sol_trap.m_f, sol_trap.fuel, t_trap);

fprintf('\n=== Nonlinear ZOH with RK4 propagation (Task 2 variant a) ===\n');
t0 = tic;
sol_zoh_nd = solve_zoh(tf_nd, N, dnd, n_sub);
t_zoh = toc(t0);
sol_zoh = dim_sol(sol_zoh_nd, ref);
fprintf('  m_f = %.2f kg, fuel = %.2f kg, wall = %.1f s\n', ...
    sol_zoh.m_f, sol_zoh.fuel, t_zoh);

fprintf('\n=== LTV-linearized ZOH with SCvx (Task 2 variant b, Appendix A) ===\n');
% Warm-start from the trapezoidal solution: the LTV linearisation is only
% locally accurate, so a near-optimal initial reference + a hard trust region
% on each decision variable are needed to keep SCvx from drifting.
t0 = tic;
[sol_scvx_nd, scvx_hist] = solve_scvx(tf_nd, N, dnd, scvx_max, scvx_tol, ...
                                      sol_trap_nd, trust);
t_scvx = toc(t0);
sol_scvx = dim_sol(sol_scvx_nd, ref);
sol_scvx.iter = sol_scvx_nd.iter;
fprintf('  m_f = %.2f kg, fuel = %.2f kg, iter = %d, wall = %.1f s\n', ...
    sol_scvx.m_f, sol_scvx.fuel, sol_scvx.iter, t_scvx);

if yalmip_ok
    fprintf('\n=== LTV-linearized ZOH with SCvx + YALMIP/ECOS (Task 2 variant c) ===\n');
    t0_y = tic;
    [sol_yscvx_nd, scvx_hist_y] = solve_scvx_yalmip(tf_nd, N, dnd, scvx_max, scvx_tol, ...
                                                     sol_trap_nd, trust);
    t_yalmip = toc(t0_y);
    sol_yscvx = dim_sol(sol_yscvx_nd, ref);
    sol_yscvx.iter = sol_yscvx_nd.iter;
    fprintf('  m_f = %.2f kg, fuel = %.2f kg, iter = %d, wall = %.2f s\n', ...
        sol_yscvx.m_f, sol_yscvx.fuel, sol_yscvx.iter, t_yalmip);
end

%% Forward-integration validation (continuous-time fidelity check)
% Performed in non-dim, error norms reported in non-dim too (state-norm).
[~, X_tfi] = fwd_integrate(sol_trap_nd, dnd, 'pwl');
[~, X_zfi] = fwd_integrate(sol_zoh_nd,  dnd, 'zoh');
[~, X_sfi] = fwd_integrate(sol_scvx_nd, dnd, 'zoh');

err_trap = node_err(sol_trap_nd, X_tfi);
err_zoh  = node_err(sol_zoh_nd,  X_zfi);
err_scvx = node_err(sol_scvx_nd, X_sfi);
if yalmip_ok
    [~, X_yfi] = fwd_integrate(sol_yscvx_nd, dnd, 'zoh');
    err_yscvx  = node_err(sol_yscvx_nd, X_yfi);
end

fprintf('\n=== Transcription fidelity (max grid-node nondim state error) ===\n');
fprintf('  Trapezoidal              : max ||delta x|| = %.4e\n', max(err_trap));
fprintf('  ZOH (RK4)                : max ||delta x|| = %.4e\n', max(err_zoh));
fprintf('  ZOH (LTV + SCvx)         : max ||delta x|| = %.4e\n', max(err_scvx));
if yalmip_ok
    fprintf('  ZOH (LTV + SCvx YALMIP)  : max ||delta x|| = %.4e\n', max(err_yscvx));
end

%% Replay landing accuracy + wall-time summary (SI units)
%  Touchdown dispersion when the optimised control schedule is flown through
%  the continuous-time nonlinear dynamics (ode45 replay): position and
%  velocity error norms at t_f w.r.t. the (0,0,0,0) target, and final-mass
%  drift w.r.t. the transcription's own m_f.
land = @(s_nd, X) [norm(X(end,1:2)) * ref.L, ...
                   norm(X(end,3:4)) * ref.V, ...
                   (s_nd.m_f - X(end,5)) * ref.m];
acc_trap = land(sol_trap_nd, X_tfi);
acc_zoh  = land(sol_zoh_nd,  X_zfi);
acc_scvx = land(sol_scvx_nd, X_sfi);

fprintf('\n=== Replay landing accuracy (ode45, SI) and wall time ===\n');
fprintf('%-26s | %12s | %14s | %14s | %9s\n', ...
    'Transcription', 'pos err [m]', 'vel err [m/s]', 'm_f drift [kg]', 'wall [s]');
fprintf('%-26s | %12.3e | %14.3e | %14.3e | %9.1f\n', ...
    'Trapezoidal (PWL)', acc_trap, t_trap);
fprintf('%-26s | %12.3e | %14.3e | %14.3e | %9.1f\n', ...
    'ZOH (RK4)', acc_zoh, t_zoh);
fprintf('%-26s | %12.3e | %14.3e | %14.3e | %9.1f\n', ...
    'ZOH (LTV + SCvx)', acc_scvx, t_scvx);
if yalmip_ok
    acc_yscvx = land(sol_yscvx_nd, X_yfi);
    fprintf('%-26s | %12.3e | %14.3e | %14.3e | %9.1f\n', ...
        'ZOH (LTV + SCvx YALMIP)', acc_yscvx, t_yalmip);
end

%% Plots
if yalmip_ok
    plot_compare4(sol_trap, sol_zoh, sol_scvx, sol_yscvx, ...
                  err_trap, err_zoh, err_scvx, err_yscvx, ...
                  scvx_hist, scvx_hist_y, data);
else
    plot_compare3(sol_trap, sol_zoh, sol_scvx, ...
                  err_trap, err_zoh, err_scvx, scvx_hist, data);
end

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
        fullfile(fig_dir, ['task2_' slugify(nm) '.png']), 'Resolution', 200);
end

%% =====================================================================
%  Local functions
%  =====================================================================

function [ref, dnd] = nondim(d)
    ref.L = d.y0;
    ref.g = d.g;
    ref.t = sqrt(ref.L / ref.g);
    ref.V = sqrt(ref.g * ref.L);
    ref.m = d.m0;
    ref.T = ref.m * ref.g;
    dnd.x0       = d.x0  / ref.L;
    dnd.y0       = d.y0  / ref.L;
    dnd.vx0      = d.vx0 / ref.V;
    dnd.vy0      = d.vy0 / ref.V;
    dnd.m0       = d.m0  / ref.m;
    dnd.Tmin     = d.Tmin / ref.T;
    dnd.Tmax     = d.Tmax / ref.T;
    dnd.Vc       = ref.V / d.c;
    dnd.theta_mx = d.theta_mx;
end

function sol = dim_sol(s_nd, ref)
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

function dx = rhs(x, u, Vc)
    % Non-dim continuous dynamics.
    Tmag = sqrt(u(1)^2 + u(2)^2);
    dx = [ x(3); x(4); u(1)/x(5); u(2)/x(5) - 1; -Vc * Tmag ];
end

function [A_jac, B_jac] = jacobians(x, u, Vc)
    % Jacobians of f w.r.t. (x, u) at point (x, u), in non-dim form.
    m = x(5);  Tx = u(1);  Ty = u(2);
    Tmag_reg = sqrt(Tx^2 + Ty^2 + 1e-6);
    A_jac = zeros(5);
    A_jac(1,3) = 1;
    A_jac(2,4) = 1;
    A_jac(3,5) = -Tx / m^2;
    A_jac(4,5) = -Ty / m^2;
    B_jac = zeros(5, 2);
    B_jac(3,1) = 1/m;
    B_jac(4,2) = 1/m;
    B_jac(5,1) = -Vc * Tx / Tmag_reg;
    B_jac(5,2) = -Vc * Ty / Tmag_reg;
end

function x_next = rk4_zoh(x, u, dt, Vc, n_sub)
    h = dt / n_sub;
    for ii = 1:n_sub
        k1 = rhs(x,                u, Vc);
        k2 = rhs(x + 0.5*h*k1,     u, Vc);
        k3 = rhs(x + 0.5*h*k2,     u, Vc);
        k4 = rhs(x +     h*k3,     u, Vc);
        x  = x + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
    end
    x_next = x;
end

function dz = ltv_aug_rhs(~, z, u_k, Vc)
    % Augmented state for the Appendix A construction.
    %   z(1:5)   = x_ref(t)
    %   z(6:30)  = vec(Phi(t,t_k))   (5x5)
    %   z(31:40) = vec(B_hat(t))     (5x2 running integral)
    %   z(41:45) = c_hat(t)          (5x1 running integral)
    x_ref = z(1:5);
    Phi   = reshape(z(6:30),  5, 5);
    Bhat  = reshape(z(31:40), 5, 2);
    chat  = z(41:45);
    [A_jac, B_jac] = jacobians(x_ref, u_k, Vc);
    f_val = rhs(x_ref, u_k, Vc);
    c_off = f_val - A_jac * x_ref - B_jac * u_k;
    dx_ref = f_val;
    dPhi   = A_jac * Phi;
    dBhat  = A_jac * Bhat + B_jac;
    dchat  = A_jac * chat + c_off;
    dz = [dx_ref; dPhi(:); dBhat(:); dchat];
end

function [Abar, Bbar, cbar] = compute_ltv_zoh(ref, tf, N, d)
    dt = tf / (N - 1);
    Abar = zeros(5, 5, N-1);
    Bbar = zeros(5, 2, N-1);
    cbar = zeros(5,    N-1);
    opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
    for k = 1:N-1
        x_k = [ref.x(k); ref.y(k); ref.vx(k); ref.vy(k); ref.m(k)];
        u_k = [ref.Tx(k); ref.Ty(k)];
        z0  = [x_k; reshape(eye(5), [], 1); zeros(10,1); zeros(5,1)];
        [~, Z] = ode45(@(t,z) ltv_aug_rhs(t, z, u_k, d.Vc), [0, dt], z0, opts);
        zf = Z(end, :).';
        Abar(:,:,k) = reshape(zf(6:30),  5, 5);
        Bbar(:,:,k) = reshape(zf(31:40), 5, 2);
        cbar(:,k)   = zf(41:45);
    end
end

function sol = solve_ltv_nlp(tf, N, d, Abar, Bbar, cbar, ref, trust)
    nz  = 7*N;
    idx = @(i) (i-1)*7 + (1:7);
    if nargin >= 7 && ~isempty(ref)
        z0 = ref_to_z(ref, idx, N);
    else
        z0 = init_guess(N, d, idx, true);
    end
    [lb, ub] = box_bounds(N, d, idx, true);
    if nargin >= 8 && ~isempty(trust) && ~isempty(ref)
        [lb, ub] = apply_trust(lb, ub, ref, idx, N, trust);
    end

    % Build linear equality block from triplets (avoids SPRIX warning)
    n_eq    = 9 + 5*(N-1);
    nnz_dyn = 8 * 5 * (N-1);
    nnz_tot = 9 + nnz_dyn;
    rows = zeros(nnz_tot, 1);
    cols = zeros(nnz_tot, 1);
    vals = zeros(nnz_tot, 1);
    beq  = zeros(n_eq, 1);

    s1 = idx(1);  sN = idx(N);
    rows(1:9) = 1:9;
    cols(1:9) = [s1(1) s1(2) s1(3) s1(4) s1(5) sN(1) sN(2) sN(3) sN(4)];
    vals(1:9) = 1;
    beq(1:5)  = [d.x0; d.y0; d.vx0; d.vy0; d.m0];

    p = 9;
    for k = 1:N-1
        sk = idx(k);  skp = idx(k+1);
        base_row = 9 + 5*(k-1);
        for i = 1:5
            row = base_row + i;
            rows(p+1) = row;  cols(p+1) = skp(i);  vals(p+1) = 1;
            rows(p+2:p+6) = row;
            cols(p+2:p+6) = sk(1:5);
            vals(p+2:p+6) = -Abar(i, 1:5, k);
            rows(p+7) = row;  cols(p+7) = sk(6);  vals(p+7) = -Bbar(i, 1, k);
            rows(p+8) = row;  cols(p+8) = sk(7);  vals(p+8) = -Bbar(i, 2, k);
            beq(row) = cbar(i, k);
            p = p + 8;
        end
    end
    Aeq = sparse(rows, cols, vals, n_eq, nz);

    iN_m = (N-1)*7 + 5;
    [z_opt, ~, ef] = fmincon(@(z) -z(iN_m), z0, [], [], full(Aeq), beq, lb, ub, ...
                             @(z) ltv_nonlcon(z, N, d), fmincon_opts('off'));
    if ef <= 0, warning('LTV fmincon ef = %d', ef); end
    sol = unpack(z_opt, tf, N, d);
end

function [c_ineq, c_eq] = ltv_nonlcon(z, N, d)
    Z = reshape(z, 7, N);
    c_ineq = path_ineq(Z, d);
    c_eq   = [];
end

function [sol, hist] = solve_scvx(tf, N, d, max_iter, tol, init_ref, base_trust)
    % SCvx outer loop with adaptive trust-region ratio.  At each iteration
    % the LTV NLP is solved within the current trust region; the candidate
    % step is then validated against the nonlinear dynamics by forward
    % integration.  The trust-region scale rho is grown/shrunk based on the
    % ratio eta = (actual cost reduction)/(predicted cost reduction):
    %   eta < eta_l       -> reject step, shrink rho
    %   eta_l <= eta < eta_h -> accept step, keep rho
    %   eta >= eta_h      -> accept step, grow rho
    % This stabilises SCvx against linearisation-error exploitation that
    % otherwise causes oscillation around the true optimum.
    idx = @(i) (i-1)*7 + (1:7);
    if nargin < 6 || isempty(init_ref)
        z0  = init_guess(N, d, idx, true);
        ref = unpack(z0, tf, N, d);
    else
        ref = init_ref;
    end
    if nargin < 7
        base_trust = struct('pos', 0.17, 'vel', 0.6, 'mass', 0.1, 'thrust', 1.0);
    end

    rho     = 1.0;       % trust-region scale (multiplies base_trust radii)
    rho_min = 1e-3;
    rho_max = 1.0;
    eta_l   = 0.25;
    eta_h   = 0.7;

    hist.m_f   = nan(max_iter, 1);
    hist.delta = nan(max_iter, 1);
    hist.rho   = nan(max_iter, 1);
    hist.eta   = nan(max_iter, 1);
    hist.acc   = false(max_iter, 1);

    sol_best = ref;
    converged = false;
    for iter = 1:max_iter
        scaled = struct('pos',    rho*base_trust.pos, ...
                        'vel',    rho*base_trust.vel, ...
                        'mass',   rho*base_trust.mass, ...
                        'thrust', rho*base_trust.thrust);

        [Abar, Bbar, cbar] = compute_ltv_zoh(ref, tf, N, d);
        sol_cand = solve_ltv_nlp(tf, N, d, Abar, Bbar, cbar, ref, scaled);

        % Trust-region ratio: predicted vs. actual cost change.  Cost is
        % -m_f, so a positive cost *reduction* corresponds to m_f going up.
        J_pred = sol_cand.m_f - ref.m_f;             % LTV-predicted gain in m_f
        [~, X_act]  = fwd_integrate(sol_cand, d, 'zoh');
        m_f_actual  = X_act(end, 5);
        J_act       = m_f_actual - ref.m_f;          % actual gain in m_f
        if abs(J_pred) < 1e-10
            eta = 1;                                 % no proposed change
        else
            eta = J_act / J_pred;
        end

        delta_x = norm( ...
            [sol_cand.x  - ref.x;  sol_cand.y  - ref.y; ...
             sol_cand.vx - ref.vx; sol_cand.vy - ref.vy; ...
             sol_cand.m  - ref.m]);

        accepted = (eta >= eta_l);
        hist.m_f(iter)   = sol_cand.m_f;
        hist.delta(iter) = delta_x;
        hist.rho(iter)   = rho;
        hist.eta(iter)   = eta;
        hist.acc(iter)   = accepted;
        fprintf('  SCvx iter %2d:  rho=%.3f  eta=%+7.3f  delta_x=%.3e  m_f=%.4f  %s\n', ...
            iter, rho, eta, delta_x, sol_cand.m_f, ...
            ternary(accepted,'ACCEPTED','rejected'));

        if accepted
            sol_best = sol_cand;
            ref = sol_cand;
            if eta > eta_h
                rho = min(rho_max, 2 * rho);
            end
            if delta_x < tol
                fprintf('    SCvx converged (delta_x < tol).\n');
                converged = true;
                break;
            end
        else
            rho = 0.5 * rho;
            if rho < rho_min
                fprintf('    SCvx stopping: trust region collapsed.\n');
                break;
            end
        end
    end
    if ~converged && iter == max_iter
        fprintf('    SCvx hit iteration cap (no further improvement).\n');
    end
    sol = sol_best;
    sol.iter = iter;
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end

function z = ref_to_z(ref, idx, N)
    z = zeros(7*N, 1);
    for i = 1:N
        s = idx(i);
        z(s(1)) = ref.x(i);
        z(s(2)) = ref.y(i);
        z(s(3)) = ref.vx(i);
        z(s(4)) = ref.vy(i);
        z(s(5)) = ref.m(i);
        z(s(6)) = ref.Tx(i);
        z(s(7)) = ref.Ty(i);
    end
end

function [lb, ub] = apply_trust(lb, ub, ref, idx, N, trust)
    % Hard box trust region centred at the reference, intersected with the
    % existing box bounds (in non-dim units).
    for i = 1:N
        s = idx(i);
        lb(s(1)) = max(lb(s(1)), ref.x(i)  - trust.pos);
        ub(s(1)) = min(ub(s(1)), ref.x(i)  + trust.pos);
        lb(s(2)) = max(lb(s(2)), ref.y(i)  - trust.pos);
        ub(s(2)) = min(ub(s(2)), ref.y(i)  + trust.pos);
        lb(s(3)) = max(lb(s(3)), ref.vx(i) - trust.vel);
        ub(s(3)) = min(ub(s(3)), ref.vx(i) + trust.vel);
        lb(s(4)) = max(lb(s(4)), ref.vy(i) - trust.vel);
        ub(s(4)) = min(ub(s(4)), ref.vy(i) + trust.vel);
        lb(s(5)) = max(lb(s(5)), ref.m(i)  - trust.mass);
        ub(s(5)) = min(ub(s(5)), ref.m(i)  + trust.mass);
        if i < N
            lb(s(6)) = max(lb(s(6)), ref.Tx(i) - trust.thrust);
            ub(s(6)) = min(ub(s(6)), ref.Tx(i) + trust.thrust);
            lb(s(7)) = max(lb(s(7)), ref.Ty(i) - trust.thrust);
            ub(s(7)) = min(ub(s(7)), ref.Ty(i) + trust.thrust);
        end
    end
end

function sol = solve_trap(tf, N, d)
    dt = tf / (N - 1);  idx = @(i) (i-1)*7 + (1:7);
    z0 = init_guess(N, d, idx, false);
    [lb, ub]   = box_bounds(N, d, idx, false);
    [Aeq, beq] = bcs(N, d, idx);
    iN_m = (N-1)*7 + 5;
    [z_opt, ~, ef] = fmincon(@(z) -z(iN_m), z0, [], [], full(Aeq), beq, lb, ub, ...
                             @(z) trap_nonlcon(z, N, dt, d), fmincon_opts());
    if ef <= 0, warning('trap fmincon ef = %d', ef); end
    sol = unpack(z_opt, tf, N, d);
end

function sol = solve_zoh(tf, N, d, n_sub)
    dt = tf / (N - 1);  idx = @(i) (i-1)*7 + (1:7);
    z0 = init_guess(N, d, idx, true);
    [lb, ub]   = box_bounds(N, d, idx, true);
    [Aeq, beq] = bcs(N, d, idx);
    iN_m = (N-1)*7 + 5;
    [z_opt, ~, ef] = fmincon(@(z) -z(iN_m), z0, [], [], full(Aeq), beq, lb, ub, ...
                             @(z) zoh_nonlcon(z, N, dt, d, n_sub), fmincon_opts());
    if ef <= 0, warning('zoh fmincon ef = %d', ef); end
    sol = unpack(z_opt, tf, N, d);
end

function z0 = init_guess(N, d, idx, zero_uN)
    z0 = zeros(7*N, 1);
    for i = 1:N
        a = (i-1) / (N-1);  s = idx(i);
        z0(s(1)) = (1-a) * d.x0;
        z0(s(2)) = (1-a) * d.y0;
        z0(s(3)) = (1-a) * d.vx0;
        z0(s(4)) = (1-a) * d.vy0;
        z0(s(5)) = d.m0 * (1 - 0.3*a);
        z0(s(6)) = 0;
        if zero_uN && i == N
            z0(s(7)) = 0;
        else
            z0(s(7)) = d.m0;     % hover (gravity = 1 in nondim)
        end
    end
end

function [lb, ub] = box_bounds(N, d, idx, zero_uN)
    nz = 7*N;  lb = -inf(nz,1);  ub = inf(nz,1);
    for i = 1:N
        s = idx(i);
        lb(s(2)) = 0;
        lb(s(5)) = 1e-3;  ub(s(5)) = d.m0;
        if zero_uN && i == N
            lb(s(6)) = 0;  ub(s(6)) = 0;
            lb(s(7)) = 0;  ub(s(7)) = 0;
        else
            lb(s(6)) = -d.Tmax;  ub(s(6)) = d.Tmax;
            lb(s(7)) = -d.Tmax;  ub(s(7)) = d.Tmax;
        end
    end
end

function [Aeq, beq] = bcs(N, d, idx)
    s1 = idx(1);  sN = idx(N);
    Aeq = sparse(9, 7*N);  beq = zeros(9,1);
    Aeq(1, s1(1)) = 1;  beq(1) = d.x0;
    Aeq(2, s1(2)) = 1;  beq(2) = d.y0;
    Aeq(3, s1(3)) = 1;  beq(3) = d.vx0;
    Aeq(4, s1(4)) = 1;  beq(4) = d.vy0;
    Aeq(5, s1(5)) = 1;  beq(5) = d.m0;
    Aeq(6, sN(1)) = 1;  Aeq(7, sN(2)) = 1;
    Aeq(8, sN(3)) = 1;  Aeq(9, sN(4)) = 1;
end

function opts = fmincon_opts(display)
    if nargin < 1, display = 'final'; end
    opts = optimoptions('fmincon', ...
        'Algorithm', 'sqp', 'Display', display, ...
        'MaxIterations', 1000, 'MaxFunctionEvaluations', 1e6, ...
        'OptimalityTolerance', 1e-5, 'ConstraintTolerance', 1e-6, ...
        'StepTolerance', 1e-10);
end

function sol = unpack(z, tf, N, d)
    Z = reshape(z, 7, N).';
    sol.t   = linspace(0, tf, N).';
    sol.x   = Z(:,1);   sol.y   = Z(:,2);
    sol.vx  = Z(:,3);   sol.vy  = Z(:,4);
    sol.m   = Z(:,5);
    sol.Tx  = Z(:,6);   sol.Ty  = Z(:,7);
    sol.Tmag = sqrt(sol.Tx.^2 + sol.Ty.^2);
    sol.tf   = tf;
    sol.m_f  = sol.m(end);
    sol.m0   = d.m0;
end

function [c_ineq, c_eq] = trap_nonlcon(z, N, dt, d)
    Z = reshape(z, 7, N);
    f = zeros(5, N);
    for i = 1:N
        f(:,i) = rhs(Z(1:5,i), Z(6:7,i), d.Vc);
    end
    defs = zeros(5, N-1);
    for k = 1:N-1
        defs(:,k) = Z(1:5, k+1) - Z(1:5, k) - 0.5*dt*(f(:,k) + f(:,k+1));
    end
    c_eq   = defs(:);
    c_ineq = path_ineq(Z, d);
end

function [c_ineq, c_eq] = zoh_nonlcon(z, N, dt, d, n_sub)
    Z = reshape(z, 7, N);
    defs = zeros(5, N-1);
    for k = 1:N-1
        x_pred = rk4_zoh(Z(1:5,k), Z(6:7,k), dt, d.Vc, n_sub);
        defs(:,k) = Z(1:5, k+1) - x_pred;
    end
    c_eq   = defs(:);
    c_ineq = path_ineq(Z, d);
end

function c_ineq = path_ineq(Z, d)
    Tmag      = sqrt(Z(6,:).^2 + Z(7,:).^2).';
    g_thr_lo  = d.Tmin - Tmag;
    g_thr_hi  = Tmag - d.Tmax;
    tt        = tan(d.theta_mx);
    g_gs_pos  = ( Z(1,:).' - tt*Z(2,:).');
    g_gs_neg  = (-Z(1,:).' - tt*Z(2,:).');
    c_ineq    = [g_thr_lo; g_thr_hi; g_gs_pos; g_gs_neg];
end

function [t, X] = fwd_integrate(sol, d, mode)
    % Forward-integrate the *continuous-time* nonlinear dynamics using the
    % optimised (non-dim) control schedule, sampling at grid nodes for
    % direct comparison with the discretised solution.
    N = numel(sol.t);
    X = zeros(N, 5);
    X(1,:) = [d.x0, d.y0, d.vx0, d.vy0, d.m0];
    opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
    for k = 1:N-1
        t_k = sol.t(k);  t_kp = sol.t(k+1);
        switch mode
            case 'zoh'
                u_fcn = @(tt) [sol.Tx(k); sol.Ty(k)];
            case 'pwl'
                u_fcn = @(tt) [
                    sol.Tx(k) + (sol.Tx(k+1)-sol.Tx(k))*(tt-t_k)/(t_kp-t_k);
                    sol.Ty(k) + (sol.Ty(k+1)-sol.Ty(k))*(tt-t_k)/(t_kp-t_k)];
            otherwise
                error('Unknown mode: %s', mode);
        end
        rhs_t = @(tt, x) rhs(x, u_fcn(tt), d.Vc);
        [~, Y] = ode45(rhs_t, [t_k, t_kp], X(k,:).', opts);
        X(k+1,:) = Y(end,:);
    end
    t = sol.t;
end

function e = node_err(sol, X)
    e = vecnorm( ...
        [sol.x sol.y sol.vx sol.vy] - X(:,1:4), 2, 2);
end

function plot_compare3(s_t, s_z, s_s, err_t, err_z, err_s, hist, d)
    cT = [0.0 0.4 0.8];
    cZ = [0.85 0.33 0.1];
    cS = [0.47 0.67 0.19];

    % --- Trajectory ---
    figure('Name','Trajectory comparison','Position',[100 100 600 500]);
    hold on; grid on; axis equal;
    yy = linspace(0, max([s_t.y; s_z.y; s_s.y])*1.05, 50);
    xx = tan(d.theta_mx) * yy;
    plot( xx, yy, '--', 'Color',[0.5 0.5 0.5], 'HandleVisibility','off');
    plot(-xx, yy, '--', 'Color',[0.5 0.5 0.5], 'HandleVisibility','off');
    plot(s_t.x, s_t.y, '-',  'Color', cT, 'LineWidth', 1.6, 'DisplayName','Trapezoidal');
    plot(s_z.x, s_z.y, '--', 'Color', cZ, 'LineWidth', 1.6, 'DisplayName','ZOH (RK4)');
    plot(s_s.x, s_s.y, ':',  'Color', cS, 'LineWidth', 2.0, 'DisplayName','ZOH (LTV+SCvx)');
    plot(0, 0, 'k^', 'MarkerSize', 8, 'MarkerFaceColor','k', 'HandleVisibility','off');
    xlabel('x  [m]'); ylabel('y  [m]');
    title('Descent trajectory: three transcriptions vs glide-slope corridor');
    legend('Location','best');

    % --- Thrust magnitude ---
    figure('Name','Thrust comparison','Position',[100 100 600 400]);
    hold on; grid on;
    plot( s_t.t, s_t.Tmag/1e3, '-',  'Color', cT, 'LineWidth', 1.6, 'DisplayName','Trapezoidal (PWL)');
    stairs(s_z.t, s_z.Tmag/1e3, '--', 'Color', cZ, 'LineWidth', 1.6, 'DisplayName','ZOH RK4 (PWC)');
    stairs(s_s.t, s_s.Tmag/1e3, ':',  'Color', cS, 'LineWidth', 2.0, 'DisplayName','ZOH LTV+SCvx (PWC)');
    yline(d.Tmax/1e3, 'k--', 'T_{max}', 'HandleVisibility','off');
    xlabel('t  [s]'); ylabel('|T|  [kN]');
    title('Thrust magnitude'); legend('Location','best');

    % --- Mass ---
    figure('Name','Mass comparison','Position',[100 100 600 400]);
    hold on; grid on;
    plot(s_t.t, s_t.m, '-',  'Color', cT, 'LineWidth', 1.6, 'DisplayName','Trapezoidal');
    plot(s_z.t, s_z.m, '--', 'Color', cZ, 'LineWidth', 1.6, 'DisplayName','ZOH (RK4)');
    plot(s_s.t, s_s.m, ':',  'Color', cS, 'LineWidth', 2.0, 'DisplayName','ZOH (LTV+SCvx)');
    xlabel('t  [s]'); ylabel('m  [kg]');
    title('Vehicle mass'); legend('Location','best');

    % --- Transcription fidelity (nondim error vs dim time) ---
    figure('Name','Transcription fidelity','Position',[100 100 600 400]);
    semilogy(s_t.t, max(err_t, eps), '-',  'Color', cT, 'LineWidth', 1.6, ...
        'DisplayName','Trapezoidal (PWL fwd-int)');
    hold on; grid on;
    semilogy(s_z.t, max(err_z, eps), '--', 'Color', cZ, 'LineWidth', 1.6, ...
        'DisplayName','ZOH RK4 (PWC fwd-int)');
    semilogy(s_s.t, max(err_s, eps), ':',  'Color', cS, 'LineWidth', 2.0, ...
        'DisplayName','ZOH LTV+SCvx (PWC fwd-int)');
    xlabel('t  [s]'); ylabel('||x_{NLP} - x_{ode45}|| (nondim)');
    title('Transcription fidelity at grid nodes (non-dim state norm)');
    legend('Location','best');

    % --- SCvx convergence trace (adaptive trust region) ---
    valid = ~isnan(hist.delta);
    iters = find(valid);
    figure('Name','SCvx convergence','Position',[100 100 800 500]);
    subplot(2,1,1);
    yyaxis left;
    semilogy(iters, max(hist.delta(valid), eps), '-o', 'LineWidth', 1.4, 'DisplayName','||\Delta x||');
    hold on;
    semilogy(iters, max(hist.rho(valid),   eps), '-s', 'LineWidth', 1.4, 'DisplayName','\rho (trust scale)');
    ylabel('log scale');
    yyaxis right;
    plot(iters, hist.m_f(valid), '-d', 'LineWidth', 1.4, 'Color',[0.4 0.4 0.4], 'DisplayName','m_f (nondim)');
    ylabel('m_f / m_0');
    grid on;
    legend('Location','best');
    title('SCvx convergence: state delta, trust scale, final mass');

    subplot(2,1,2);
    bar_h = bar(iters, hist.eta(valid), 0.7);
    ylim([-0.5 2]);
    yline(0.25, 'k--', '\eta_l = 0.25');
    yline(0.7,  'k--', '\eta_h = 0.7');
    if isfield(hist, 'acc')
        for k = 1:length(iters)
            if ~hist.acc(iters(k))
                bar_h.FaceColor = 'flat';
                bar_h.CData(k,:) = [0.85 0.33 0.1];
            end
        end
    end
    xlabel('SCvx iteration');
    ylabel('Trust ratio \eta');
    grid on;
    title('Trust-region ratio \eta = (actual gain)/(predicted gain) per iteration');
end

function sol = solve_ltv_nlp_yalmip(tf, N, d, Abar, Bbar, cbar, ref_sol, trust)
% Inner SOCP for one SCvx iteration, solved via YALMIP + ECOS.
% The LTV dynamics are linear equalities; thrust upper-bound is a
% second-order cone constraint; glide-slope and trust-region are linear.
    X = sdpvar(5, N,   'full');   % state  [x; y; vx; vy; m] at each node
    U = sdpvar(2, N-1, 'full');   % control [Tx; Ty] per ZOH interval

    tt = tan(d.theta_mx);

    cstr = [X(:,1) == [d.x0; d.y0; d.vx0; d.vy0; d.m0]];   % I.C.
    cstr = [cstr, X(1:4, N) == 0];                            % terminal BCs

    for k = 1:N-1   % LTV dynamics (linear equalities)
        cstr = [cstr, X(:,k+1) == Abar(:,:,k)*X(:,k) + Bbar(:,:,k)*U(:,k) + cbar(:,k)];
    end

    for k = 1:N-1   % thrust upper bound (SOCP); Tmin=0 is trivially satisfied
        cstr = [cstr, norm(U(:,k)) <= d.Tmax];
    end

    for k = 1:N     % glide-slope cone (linear)
        cstr = [cstr,  X(1,k) <= tt * X(2,k)];
        cstr = [cstr, -X(1,k) <= tt * X(2,k)];
    end

    cstr = [cstr, X(2,:) >= 0];                          % altitude >= 0
    cstr = [cstr, X(5,:) >= 1e-3, X(5,:) <= d.m0];      % mass bounds

    if nargin >= 7 && ~isempty(ref_sol) && nargin >= 8 && ~isempty(trust)
        for k = 1:N
            cstr = [cstr, X(1,k) >= ref_sol.x(k)  - trust.pos,  X(1,k) <= ref_sol.x(k)  + trust.pos];
            cstr = [cstr, X(2,k) >= ref_sol.y(k)  - trust.pos,  X(2,k) <= ref_sol.y(k)  + trust.pos];
            cstr = [cstr, X(3,k) >= ref_sol.vx(k) - trust.vel,  X(3,k) <= ref_sol.vx(k) + trust.vel];
            cstr = [cstr, X(4,k) >= ref_sol.vy(k) - trust.vel,  X(4,k) <= ref_sol.vy(k) + trust.vel];
            cstr = [cstr, X(5,k) >= ref_sol.m(k)  - trust.mass, X(5,k) <= ref_sol.m(k)  + trust.mass];
        end
        for k = 1:N-1
            cstr = [cstr, U(1,k) >= ref_sol.Tx(k) - trust.thrust, U(1,k) <= ref_sol.Tx(k) + trust.thrust];
            cstr = [cstr, U(2,k) >= ref_sol.Ty(k) - trust.thrust, U(2,k) <= ref_sol.Ty(k) + trust.thrust];
        end
    end

    res = optimize(cstr, -X(5,N), sdpsettings('solver','ecos','verbose',0));
    if res.problem ~= 0
        warning('YALMIP/ECOS: problem flag = %d (%s)', res.problem, res.info);
    end

    X_val = value(X);
    U_val = [value(U), zeros(2,1)];   % pad: last node has no ZOH control
    sol.t    = linspace(0, tf, N).';
    sol.x    = X_val(1,:).';   sol.y  = X_val(2,:).';
    sol.vx   = X_val(3,:).';   sol.vy = X_val(4,:).';
    sol.m    = X_val(5,:).';
    sol.Tx   = U_val(1,:).';   sol.Ty = U_val(2,:).';
    sol.Tmag = sqrt(sol.Tx.^2 + sol.Ty.^2);
    sol.tf   = tf;   sol.m_f = sol.m(end);   sol.m0 = d.m0;
end

function [sol, hist] = solve_scvx_yalmip(tf, N, d, max_iter, tol, init_ref, base_trust)
% SCvx outer loop — same adaptive trust-region logic as solve_scvx;
% inner subproblem solved as an SOCP via solve_ltv_nlp_yalmip instead of fmincon.
    idx = @(i) (i-1)*7 + (1:7);
    if nargin < 6 || isempty(init_ref)
        ref = unpack(init_guess(N, d, idx, true), tf, N, d);
    else
        ref = init_ref;
    end
    if nargin < 7
        base_trust = struct('pos', 0.17, 'vel', 0.6, 'mass', 0.1, 'thrust', 1.0);
    end

    rho = 1.0;   rho_min = 1e-3;   rho_max = 1.0;
    eta_l = 0.25;   eta_h = 0.7;

    hist.m_f   = nan(max_iter, 1);   hist.delta = nan(max_iter, 1);
    hist.rho   = nan(max_iter, 1);   hist.eta   = nan(max_iter, 1);
    hist.acc   = false(max_iter, 1);

    sol_best  = ref;   converged = false;
    for iter = 1:max_iter
        scaled = struct('pos',    rho*base_trust.pos,    'vel',    rho*base_trust.vel, ...
                        'mass',   rho*base_trust.mass,   'thrust', rho*base_trust.thrust);

        [Abar, Bbar, cbar] = compute_ltv_zoh(ref, tf, N, d);
        sol_cand = solve_ltv_nlp_yalmip(tf, N, d, Abar, Bbar, cbar, ref, scaled);

        J_pred = sol_cand.m_f - ref.m_f;
        [~, X_act] = fwd_integrate(sol_cand, d, 'zoh');
        J_act = X_act(end, 5) - ref.m_f;
        eta   = ternary(abs(J_pred) < 1e-10, 1, J_act / J_pred);

        delta_x = norm([sol_cand.x - ref.x; sol_cand.y - ref.y; ...
                        sol_cand.vx - ref.vx; sol_cand.vy - ref.vy; sol_cand.m - ref.m]);

        accepted = (eta >= eta_l);
        hist.m_f(iter) = sol_cand.m_f;   hist.delta(iter) = delta_x;
        hist.rho(iter) = rho;             hist.eta(iter)   = eta;
        hist.acc(iter) = accepted;
        fprintf('  SCvx-Y iter %2d:  rho=%.3f  eta=%+7.3f  delta_x=%.3e  m_f=%.4f  %s\n', ...
            iter, rho, eta, delta_x, sol_cand.m_f, ternary(accepted,'ACCEPTED','rejected'));

        if accepted
            sol_best = sol_cand;   ref = sol_cand;
            if eta > eta_h,  rho = min(rho_max, 2*rho);  end
            if delta_x < tol
                fprintf('    SCvx-Y converged (delta_x < tol).\n');
                converged = true;   break;
            end
        else
            rho = 0.5 * rho;
            if rho < rho_min
                fprintf('    SCvx-Y stopping: trust region collapsed.\n');   break;
            end
        end
    end
    if ~converged && iter == max_iter
        fprintf('    SCvx-Y hit iteration cap (no further improvement).\n');
    end
    sol = sol_best;   sol.iter = iter;
end

function plot_compare4(s_t, s_z, s_s, s_y, err_t, err_z, err_s, err_y, hist_s, hist_y, d)
    cT = [0.0  0.4  0.8];
    cZ = [0.85 0.33 0.1];
    cS = [0.47 0.67 0.19];
    cY = [0.49 0.18 0.56];   % purple — YALMIP variant

    % --- Trajectory ---
    figure('Name','Trajectory comparison','Position',[100 100 600 500]);
    hold on; grid on; axis equal;
    yy = linspace(0, max([s_t.y; s_z.y; s_s.y; s_y.y])*1.05, 50);
    xx = tan(d.theta_mx) * yy;
    plot( xx, yy, '--', 'Color',[0.5 0.5 0.5], 'HandleVisibility','off');
    plot(-xx, yy, '--', 'Color',[0.5 0.5 0.5], 'HandleVisibility','off');
    plot(s_t.x, s_t.y, '-',   'Color', cT, 'LineWidth', 1.6, 'DisplayName','Trapezoidal');
    plot(s_z.x, s_z.y, '--',  'Color', cZ, 'LineWidth', 1.6, 'DisplayName','ZOH (RK4)');
    plot(s_s.x, s_s.y, ':',   'Color', cS, 'LineWidth', 2.0, 'DisplayName','ZOH (LTV+SCvx)');
    plot(s_y.x, s_y.y, '-.',  'Color', cY, 'LineWidth', 2.0, 'DisplayName','ZOH (LTV+SCvx YALMIP)');
    plot(0, 0, 'k^', 'MarkerSize', 8, 'MarkerFaceColor','k', 'HandleVisibility','off');
    xlabel('x  [m]'); ylabel('y  [m]');
    title('Descent trajectory: four transcriptions vs glide-slope corridor');
    legend('Location','best');

    % --- Thrust magnitude ---
    figure('Name','Thrust comparison','Position',[100 100 600 400]);
    hold on; grid on;
    plot(  s_t.t, s_t.Tmag/1e3, '-',  'Color', cT, 'LineWidth', 1.6, 'DisplayName','Trapezoidal (PWL)');
    stairs(s_z.t, s_z.Tmag/1e3, '--', 'Color', cZ, 'LineWidth', 1.6, 'DisplayName','ZOH RK4 (PWC)');
    stairs(s_s.t, s_s.Tmag/1e3, ':',  'Color', cS, 'LineWidth', 2.0, 'DisplayName','ZOH LTV+SCvx (PWC)');
    stairs(s_y.t, s_y.Tmag/1e3, '-.', 'Color', cY, 'LineWidth', 2.0, 'DisplayName','ZOH LTV+SCvx YALMIP (PWC)');
    yline(d.Tmax/1e3, 'k--', 'T_{max}', 'HandleVisibility','off');
    xlabel('t  [s]'); ylabel('|T|  [kN]');
    title('Thrust magnitude'); legend('Location','best');

    % --- Mass ---
    figure('Name','Mass comparison','Position',[100 100 600 400]);
    hold on; grid on;
    plot(s_t.t, s_t.m, '-',  'Color', cT, 'LineWidth', 1.6, 'DisplayName','Trapezoidal');
    plot(s_z.t, s_z.m, '--', 'Color', cZ, 'LineWidth', 1.6, 'DisplayName','ZOH (RK4)');
    plot(s_s.t, s_s.m, ':',  'Color', cS, 'LineWidth', 2.0, 'DisplayName','ZOH (LTV+SCvx)');
    plot(s_y.t, s_y.m, '-.', 'Color', cY, 'LineWidth', 2.0, 'DisplayName','ZOH (LTV+SCvx YALMIP)');
    xlabel('t  [s]'); ylabel('m  [kg]');
    title('Vehicle mass'); legend('Location','best');

    % --- Transcription fidelity ---
    figure('Name','Transcription fidelity','Position',[100 100 600 400]);
    semilogy(s_t.t, max(err_t, eps), '-',  'Color', cT, 'LineWidth', 1.6, ...
        'DisplayName','Trapezoidal (PWL fwd-int)');
    hold on; grid on;
    semilogy(s_z.t, max(err_z, eps), '--', 'Color', cZ, 'LineWidth', 1.6, ...
        'DisplayName','ZOH RK4 (PWC fwd-int)');
    semilogy(s_s.t, max(err_s, eps), ':',  'Color', cS, 'LineWidth', 2.0, ...
        'DisplayName','ZOH LTV+SCvx (PWC fwd-int)');
    semilogy(s_y.t, max(err_y, eps), '-.', 'Color', cY, 'LineWidth', 2.0, ...
        'DisplayName','ZOH LTV+SCvx YALMIP (PWC fwd-int)');
    xlabel('t  [s]'); ylabel('||x_{NLP} - x_{ode45}|| (nondim)');
    title('Transcription fidelity at grid nodes (non-dim state norm)');
    legend('Location','best');

    % --- SCvx convergence: fmincon vs YALMIP side by side ---
    for ii = 1:2
        if ii == 1,  hist = hist_s;  lbl = 'fmincon/SQP';  col = cS;
        else,        hist = hist_y;  lbl = 'YALMIP/ECOS';  col = cY;
        end
        valid = ~isnan(hist.delta);
        iters = find(valid);
        figure('Name', sprintf('SCvx convergence — %s', lbl), 'Position',[100+400*(ii-1) 600 800 500]);
        subplot(2,1,1);
        yyaxis left;
        semilogy(iters, max(hist.delta(valid), eps), '-o', 'LineWidth', 1.4, ...
            'Color', col, 'DisplayName','||\Delta x||');
        hold on;
        semilogy(iters, max(hist.rho(valid), eps), '-s', 'LineWidth', 1.4, ...
            'Color', col*0.6, 'DisplayName','\rho (trust scale)');
        ylabel('log scale');
        yyaxis right;
        plot(iters, hist.m_f(valid), '-d', 'LineWidth', 1.4, 'Color',[0.4 0.4 0.4], ...
            'DisplayName','m_f (nondim)');
        ylabel('m_f / m_0');
        grid on; legend('Location','best');
        title(sprintf('SCvx convergence [%s]: state delta, trust scale, final mass', lbl));
        subplot(2,1,2);
        bar_h = bar(iters, hist.eta(valid), 0.7);
        bar_h.FaceColor = col;
        ylim([-0.5 2]);
        yline(0.25, 'k--', '\eta_l = 0.25');
        yline(0.7,  'k--', '\eta_h = 0.7');
        if isfield(hist, 'acc')
            for k = 1:length(iters)
                if ~hist.acc(iters(k))
                    bar_h.FaceColor = 'flat';
                    bar_h.CData(k,:) = [0.85 0.33 0.1];
                end
            end
        end
        xlabel('SCvx iteration');
        ylabel('Trust ratio \eta');
        grid on;
        title(sprintf('Trust-region ratio \\eta per iteration [%s]', lbl));
    end
end
