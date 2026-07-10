%% HM2 - Task 2: ZOH discretization (Optional)
%  Four ZOH transcriptions:
%    (a) Nonlinear ZOH + RK4 -- multiple-shooting NLP, x_{k+1} = RK4(x_k, u_k, dt)
%    (b) LTV ZOH (Appendix A) + SCvx outer loop -- linearize about a
%        reference, build discrete matrices from the Appendix A ODE, solve
%        the LTV NLP, update reference, repeat.
%    (c) same LTV+SCvx loop, inner step modelled as an SOCP (YALMIP + ECOS).
%    (d) GFOLD log-mass change of variables (z=ln m, u=T/m, slack sigma):
%        the dynamics become exactly LTI so the appendix-A ZOH is a single
%        matrix exponential; only the thrust upper bound is linearised and
%        iterated in the SCvx loop (SOCP via YALMIP + ECOS).
%
%  Variants (a)-(c) warm-start from the Task 1 trapezoidal baseline. Solved non-dim
%  (L_ref = y0, a_ref = g, t_ref = sqrt(L_ref/g), V_ref = sqrt(g*L_ref),
%  m_ref = m0, T_ref = m0*g; V_c = V_ref/c). Scaled back to SI for output.
%
%  Validation: replay the ZOH controls through ode45 vs the discretized
%  trajectory.
%
%  Reference: Homework 2 - Powered Descent Landing (Zavoli, April 2026)
%  Solver: fmincon (sqp). No external dependency.

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
scvx_max = 15;     % SCvx max iterations (variant b)
scvx_tol = 1e-3;   % SCvx tolerance, nondim (variant b)

%% Non-dimensionalisation
[ref, dnd] = nondim(data);
fprintf('Non-dim reference scales:\n');
fprintf('  L_ref = %.1f m,  V_ref = %.2f m/s,  t_ref = %.3f s\n', ref.L, ref.V, ref.t);
fprintf('  m_ref = %.0f kg, T_ref = %.0f N,  V_c = V_ref/c = %.4f\n', ref.m, ref.T, dnd.Vc);
tf_nd = tf / ref.t;

% SCvx trust-region radii (nondim)
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
% locally accurate, so a near-optimal reference plus a hard per-variable
% trust region keep SCvx from drifting.
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

    fprintf('\n=== GFOLD log-mass change of variables + SCvx (Task 2 variant d) ===\n');
    % z = ln(m), u = T/m, slack sigma make the dynamics exactly LTI: the
    % appendix-A ZOH is one matrix exponential and only the thrust upper
    % bound (||u|| <= Tmax*e^{-z}) is linearised and iterated. Self-starting
    % from the analytic max-thrust mass profile -- no trapezoidal warm start.
    t0_g = tic;
    [sol_gfold_nd, scvx_hist_g] = solve_gfold_scvx(tf_nd, N, dnd, scvx_max, scvx_tol);
    t_gfold = toc(t0_g);
    sol_gfold = dim_sol(sol_gfold_nd, ref);
    sol_gfold.iter = sol_gfold_nd.iter;
    fprintf('  m_f = %.2f kg, fuel = %.2f kg, iter = %d, wall = %.2f s\n', ...
        sol_gfold.m_f, sol_gfold.fuel, sol_gfold.iter, t_gfold);
end

%% Forward-integration validation (continuous-time fidelity check)
% Non-dim throughout; error norms reported as nondim state-norm.
[~, X_tfi] = fwd_integrate(sol_trap_nd, dnd, 'pwl');
[~, X_zfi] = fwd_integrate(sol_zoh_nd,  dnd, 'zoh');
[~, X_sfi] = fwd_integrate(sol_scvx_nd, dnd, 'zoh');

err_trap = node_err(sol_trap_nd, X_tfi);
err_zoh  = node_err(sol_zoh_nd,  X_zfi);
err_scvx = node_err(sol_scvx_nd, X_sfi);
if yalmip_ok
    [~, X_yfi] = fwd_integrate(sol_yscvx_nd, dnd, 'zoh');
    err_yscvx  = node_err(sol_yscvx_nd, X_yfi);
    % GFOLD holds the acceleration u (not the thrust T) piecewise constant,
    % so it is replayed with the u-ZOH convention (T = m(t)*u floats).
    [~, X_gfi] = fwd_integrate_uacc(sol_gfold_nd, dnd);
    err_gfold  = node_err(sol_gfold_nd, X_gfi);
end

fprintf('\n=== Transcription fidelity (max grid-node nondim state error) ===\n');
fprintf('  Trapezoidal              : max ||delta x|| = %.4e\n', max(err_trap));
fprintf('  ZOH (RK4)                : max ||delta x|| = %.4e\n', max(err_zoh));
fprintf('  ZOH (LTV + SCvx)         : max ||delta x|| = %.4e\n', max(err_scvx));
if yalmip_ok
    fprintf('  ZOH (LTV + SCvx YALMIP)  : max ||delta x|| = %.4e\n', max(err_yscvx));
    fprintf('  ZOH (GFOLD log-mass)     : max ||delta x|| = %.4e\n', max(err_gfold));
end

%% Replay landing accuracy + wall-time summary (SI units)
%  Touchdown dispersion from the ode45 replay: pos/vel error norms at t_f
%  vs the (0,0,0,0) target, plus final-mass drift vs the transcription's m_f.
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
    acc_gfold = land(sol_gfold_nd, X_gfi);
    fprintf('%-26s | %12.3e | %14.3e | %14.3e | %9.1f\n', ...
        'ZOH (GFOLD log-mass)', acc_gfold, t_gfold);
end

%% Plots
if yalmip_ok
    plot_compare5(sol_trap, sol_zoh, sol_scvx, sol_yscvx, sol_gfold, ...
                  err_trap, err_zoh, err_scvx, err_yscvx, err_gfold, ...
                  scvx_hist, scvx_hist_y, scvx_hist_g, data);
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
    nm = fig_handles(kk).Name;
    if isempty(nm); nm = sprintf('fig%d', kk); end
    try
        theme(fig_handles(kk), 'light');    % force light theme (ignore desktop dark mode)
        drawnow;
    catch
        fig_handles(kk).Color = 'w';        % fallback for pre-R2025a MATLAB
    end
    exportgraphics(fig_handles(kk), ...
        fullfile(fig_dir, ['task2_' slugify(nm) '.png']), 'Resolution', 200);
end

%% =====================================================================
%  Local functions
%  Continuous dynamics and ZOH propagator live in ode_descent.m / rk4_zoh.m,
%  shared with main_task1.m and the test suite.
%  Hot-loop functions (jacobians, ltv_aug_rhs, *_nonlcon, path_ineq) skip
%  arguments validation -- they sit inside the fmincon/ode45 inner loop.
%  =====================================================================

function [ref, dnd] = nondim(d)
    % Build reference scales and non-dim problem data.
    %   INPUT
    %     d   - struct: SI problem data
    %   OUTPUT
    %     ref - struct: reference scales
    %     dnd - struct: non-dim data
    arguments
        d (1,1) struct
    end
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
    % Scale a non-dim solution struct back to SI.
    %   INPUT
    %     s_nd - struct: non-dim solution
    %     ref  - struct: reference scales
    %   OUTPUT
    %     sol  - struct: SI solution
    arguments
        s_nd (1,1) struct
        ref  (1,1) struct
    end
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

function [A_jac, B_jac] = jacobians(x, u, Vc)
    % Jacobians df/dx, df/du at (x, u), non-dim.
    %   INPUT
    %     x     - state [x; y; vx; vy; m]
    %     u     - thrust [Tx; Ty]
    %     Vc    - V_ref/c
    %   OUTPUT
    %     A_jac - df/dx (5x5)
    %     B_jac - df/du (5x2)
    % No arguments validation by design: SCvx/ode45 hot loop.
    m = x(5);  Tx = u(1);  Ty = u(2);
    Tmag_reg = sqrt(Tx^2 + Ty^2 + 1e-6);   % regularize |T| to keep df/du finite at T=0
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

function dz = ltv_aug_rhs(~, z, u_k, Vc)
    % Appendix A augmented RHS, beta-gamma form: state + STM Phi + inverse STM
    % Psi + ZOH integrals Beta, Gamma (all referenced to the interval start t_k).
    %   INPUT
    %     z   - [x_ref(5); vec Phi(25); vec Psi(25); vec Beta(10); Gamma(5)]
    %     u_k - ZOH thrust over the interval
    %     Vc  - V_ref/c
    %   OUTPUT
    %     dz  - augmented derivative (70x1)
    % No arguments validation by design: SCvx/ode45 hot loop.
    x_ref = z(1:5);
    Phi   = reshape(z(6:30),  5, 5);
    Psi   = reshape(z(31:55), 5, 5);   % Psi = Phi^{-1}, inverse transition referenced to t_k
    % Beta = z(56:65), Gamma = z(66:70) are pure integrals: their derivatives
    % dBeta = Psi*B, dGamma = Psi*c do not depend on Beta, Gamma themselves.
    [A_jac, B_jac] = jacobians(x_ref, u_k, Vc);
    f_val = ode_descent(x_ref, u_k, Vc);
    c_off = f_val - A_jac * x_ref - B_jac * u_k;   % affine offset of the linearization
    dx_ref =  f_val;
    dPhi   =  A_jac * Phi;
    dPsi   = -Psi * A_jac;             % d(Phi^{-1})/dt = -Phi^{-1} A
    dBeta  =  Psi * B_jac;
    dGamma =  Psi * c_off;
    dz = [dx_ref; dPhi(:); dPsi(:); dBeta(:); dGamma];
end

function [Abar, Bbar, cbar] = compute_ltv_zoh(ref, tf, N, d)
    % Discrete LTV ZOH matrices per interval, integrating the Appendix A ODE.
    %   INPUT
    %     ref  - struct: reference trajectory
    %     tf   - flight time (nondim)
    %     N    - node count
    %     d    - struct: non-dim data
    %   OUTPUT
    %     Abar - state transition (5x5x(N-1))
    %     Bbar - input matrix (5x2x(N-1))
    %     cbar - affine term (5x(N-1))
    arguments
        ref (1,1) struct
        tf  (1,1) double {mustBePositive, mustBeFinite}
        N   (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(N, 2)}
        d   (1,1) struct
    end
    dt = tf / (N - 1);
    Abar = zeros(5, 5, N-1);
    Bbar = zeros(5, 2, N-1);
    cbar = zeros(5,    N-1);
    opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
    for k = 1:N-1
        x_k = [ref.x(k); ref.y(k); ref.vx(k); ref.vy(k); ref.m(k)];
        u_k = [ref.Tx(k); ref.Ty(k)];
        % Phi and Psi=Phi^{-1} both start at I; Beta, Gamma start at 0.
        z0  = [x_k; reshape(eye(5), [], 1); reshape(eye(5), [], 1); zeros(10,1); zeros(5,1)];
        [~, Z] = ode45(@(t,z) ltv_aug_rhs(t, z, u_k, d.Vc), [0, dt], z0, opts);
        zf = Z(end, :).';
        Phi_f = reshape(zf(6:30), 5, 5);
        Abar(:,:,k) = Phi_f;
        Bbar(:,:,k) = Phi_f * reshape(zf(56:65), 5, 2);   % B_k = Phi(t_{k+1}) * Beta(t_{k+1})
        cbar(:,k)   = Phi_f * zf(66:70);                  % c_k = Phi(t_{k+1}) * Gamma(t_{k+1})
    end
end

function sol = solve_ltv_nlp(tf, N, d, Abar, Bbar, cbar, ref, trust)
    % One SCvx inner NLP: LTV dynamics as linear equalities, fmincon.
    %   INPUT
    %     tf         - flight time (nondim)
    %     N          - node count
    %     d          - struct: non-dim data
    %     Abar,Bbar,cbar - discrete LTV matrices
    %     ref        - struct: reference (trust-region centre), optional
    %     trust      - struct: trust-region radii, optional
    %   OUTPUT
    %     sol        - struct: non-dim solution
    arguments
        tf   (1,1) double {mustBePositive, mustBeFinite}
        N    (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(N, 2)}
        d    (1,1) struct
        Abar (5,5,:) double
        Bbar (5,2,:) double
        cbar (5,:)   double
        ref   = []
        trust = []
    end
    nz  = 7*N;
    idx = @(i) (i-1)*7 + (1:7);
    if ~isempty(ref)
        z0 = ref_to_z(ref, idx, N);
    else
        z0 = init_guess(N, d, idx, true);
    end
    [lb, ub] = box_bounds(N, d, idx, true);
    if ~isempty(trust) && ~isempty(ref)
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
    % Path constraints only; LTV dynamics handled as linear equalities.
    %   INPUT
    %     z      - decision vector (7N x1)
    %     N      - node count
    %     d      - struct: non-dim data
    %   OUTPUT
    %     c_ineq - path constraints (<=0)
    %     c_eq   - empty
    % No arguments validation by design: fmincon nonlcon hot loop.
    Z = reshape(z, 7, N);
    c_ineq = path_ineq(Z, d);
    c_eq   = [];
end

function [sol, conv_hist] = solve_scvx(tf, N, d, max_iter, tol, init_ref, base_trust)
    % SCvx outer loop, fmincon inner NLP, adaptive trust region.
    %   INPUT
    %     tf         - flight time (nondim)
    %     N          - node count
    %     d          - struct: non-dim data
    %     max_iter   - iteration cap
    %     tol        - convergence tol on ||delta x||
    %     init_ref   - struct: initial reference, optional
    %     base_trust - struct: base trust radii
    %   OUTPUT
    %     sol        - struct: best non-dim solution (+ iter count)
    %     conv_hist  - struct: per-iteration trace
    % Each iter solves the LTV NLP in the trust region, then validates the
    % step against the nonlinear dynamics by forward integration. rho is
    % resized from eta = (actual gain)/(predicted gain):
    %   eta < eta_l            -> reject, shrink rho
    %   eta_l <= eta < eta_h   -> accept, keep rho
    %   eta >= eta_h           -> accept, grow rho
    arguments
        tf       (1,1) double {mustBePositive, mustBeFinite}
        N        (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(N, 2)}
        d        (1,1) struct
        max_iter (1,1) double {mustBeInteger, mustBePositive}
        tol      (1,1) double {mustBePositive}
        init_ref         = []
        base_trust (1,1) struct = struct('pos', 0.17, 'vel', 0.6, 'mass', 0.1, 'thrust', 1.0)
    end
    idx = @(i) (i-1)*7 + (1:7);
    if isempty(init_ref)
        ref = unpack(init_guess(N, d, idx, true), tf, N, d);
    else
        ref = init_ref;
    end

    rho     = 1.0;       % trust-region scale (multiplies base_trust radii)
    rho_min = 1e-3;
    rho_max = 1.0;
    eta_l   = 0.25;
    eta_h   = 0.7;

    conv_hist.m_f   = nan(max_iter, 1);
    conv_hist.delta = nan(max_iter, 1);
    conv_hist.rho   = nan(max_iter, 1);
    conv_hist.eta   = nan(max_iter, 1);
    conv_hist.acc   = false(max_iter, 1);

    sol_best = ref;
    converged = false;
    for iter = 1:max_iter
        scaled = struct('pos',    rho*base_trust.pos, ...
                        'vel',    rho*base_trust.vel, ...
                        'mass',   rho*base_trust.mass, ...
                        'thrust', rho*base_trust.thrust);

        [Abar, Bbar, cbar] = compute_ltv_zoh(ref, tf, N, d);
        sol_cand = solve_ltv_nlp(tf, N, d, Abar, Bbar, cbar, ref, scaled);

        % Trust ratio. Cost is -m_f, so a cost reduction means m_f up.
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
        conv_hist.m_f(iter)   = sol_cand.m_f;
        conv_hist.delta(iter) = delta_x;
        conv_hist.rho(iter)   = rho;
        conv_hist.eta(iter)   = eta;
        conv_hist.acc(iter)   = accepted;
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
    % Inline if-else.
    %   INPUT
    %     cond - logical
    %     a    - value if true
    %     b    - value if false
    %   OUTPUT
    %     out  - a or b
    if cond, out = a; else, out = b; end
end

function z = ref_to_z(ref, idx, N)
    % Pack a reference struct into the stacked decision vector.
    %   INPUT
    %     ref - struct: trajectory
    %     idx - node-index handle
    %     N   - node count
    %   OUTPUT
    %     z   - decision vector (7N x1)
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
    % Intersect box bounds with a trust region centred at the reference.
    %   INPUT
    %     lb,ub - existing bounds
    %     ref   - struct: trust-region centre
    %     idx   - node-index handle
    %     N     - node count
    %     trust - struct: radii (pos,vel,mass,thrust)
    %   OUTPUT
    %     lb,ub - tightened bounds
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
    % Trapezoidal collocation NLP (Task 1 baseline).
    %   INPUT
    %     tf  - flight time (nondim)
    %     N   - node count
    %     d   - struct: non-dim data
    %   OUTPUT
    %     sol - struct: non-dim solution
    arguments
        tf (1,1) double {mustBePositive, mustBeFinite}
        N  (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(N, 2)}
        d  (1,1) struct
    end
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
    % Nonlinear ZOH multiple-shooting NLP (RK4 defects).
    %   INPUT
    %     tf    - flight time (nondim)
    %     N     - node count
    %     d     - struct: non-dim data
    %     n_sub - RK4 substeps per interval
    %   OUTPUT
    %     sol   - struct: non-dim solution
    arguments
        tf    (1,1) double {mustBePositive, mustBeFinite}
        N     (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(N, 2)}
        d     (1,1) struct
        n_sub (1,1) double {mustBeInteger, mustBePositive}
    end
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
    % Linear state interp, hover thrust; zero control at last node if zero_uN.
    %   INPUT
    %     N       - node count
    %     d       - struct: non-dim data
    %     idx     - node-index handle
    %     zero_uN - force u(N)=0 (ZOH variants)
    %   OUTPUT
    %     z0      - decision vector (7N x1)
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
            z0(s(7)) = d.m0;     % hover (gravity = 1 nondim)
        end
    end
end

function [lb, ub] = box_bounds(N, d, idx, zero_uN)
    % Box bounds: y>=0, mass in [1e-3, m0], |T| components in [-Tmax, Tmax].
    %   INPUT
    %     N       - node count
    %     d       - struct: non-dim data
    %     idx     - node-index handle
    %     zero_uN - pin control to 0 at last node (ZOH variants)
    %   OUTPUT
    %     lb,ub   - bounds (7N x1)
    nz = 7*N;  lb = -inf(nz,1);  ub = inf(nz,1);
    for i = 1:N
        s = idx(i);
        lb(s(2)) = 0;
        lb(s(5)) = 1e-3;  ub(s(5)) = d.m0;   % mass strictly positive
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
    % Linear boundary conditions: full state at node 1, pos+vel zero at node N.
    %   INPUT
    %     N   - node count
    %     d   - struct: non-dim data
    %     idx - node-index handle
    %   OUTPUT
    %     Aeq - equality matrix (9 x 7N, sparse)
    %     beq - equality RHS (9x1)
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

function opts = fmincon_opts(display_mode)
    % Shared SQP fmincon options.
    %   INPUT
    %     display_mode - Display value (default 'final')
    %   OUTPUT
    %     opts         - optimoptions
    arguments
        display_mode {mustBeTextScalar} = 'final'
    end
    opts = optimoptions('fmincon', ...
        'Algorithm', 'sqp', 'Display', display_mode, ...
        'MaxIterations', 1000, 'MaxFunctionEvaluations', 1e6, ...
        'OptimalityTolerance', 1e-5, 'ConstraintTolerance', 1e-6, ...
        'StepTolerance', 1e-10);
end

function sol = unpack(z, tf, N, d)
    % Unpack a decision vector into a non-dim solution struct.
    %   INPUT
    %     z   - decision vector (7N x1)
    %     tf  - flight time (nondim)
    %     N   - node count
    %     d   - struct: non-dim data
    %   OUTPUT
    %     sol - struct: non-dim solution
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
    % Trapezoidal defects + path constraints.
    %   INPUT
    %     z      - decision vector (7N x1)
    %     N      - node count
    %     dt     - node spacing (nondim)
    %     d      - struct: non-dim data
    %   OUTPUT
    %     c_ineq - path constraints (<=0)
    %     c_eq   - collocation defects (=0)
    % No arguments validation by design: fmincon nonlcon hot loop.
    Z = reshape(z, 7, N);
    f = zeros(5, N);
    for i = 1:N
        f(:,i) = ode_descent(Z(1:5,i), Z(6:7,i), d.Vc);
    end
    defs = zeros(5, N-1);
    for k = 1:N-1
        defs(:,k) = Z(1:5, k+1) - Z(1:5, k) - 0.5*dt*(f(:,k) + f(:,k+1));
    end
    c_eq   = defs(:);
    c_ineq = path_ineq(Z, d);
end

function [c_ineq, c_eq] = zoh_nonlcon(z, N, dt, d, n_sub)
    % ZOH shooting defects (x_{k+1} - RK4) + path constraints.
    %   INPUT
    %     z      - decision vector (7N x1)
    %     N      - node count
    %     dt     - node spacing (nondim)
    %     d      - struct: non-dim data
    %     n_sub  - RK4 substeps
    %   OUTPUT
    %     c_ineq - path constraints (<=0)
    %     c_eq   - shooting defects (=0)
    % No arguments validation by design: fmincon nonlcon hot loop.
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
    % Thrust magnitude bounds + glide-slope cone (stacked).
    %   INPUT
    %     Z      - state/control matrix (7xN)
    %     d      - struct: non-dim data
    %   OUTPUT
    %     c_ineq - [thr_lo; thr_hi; gs_pos; gs_neg] (<=0)
    % No arguments validation by design: fmincon nonlcon hot loop.
    Tmag      = sqrt(Z(6,:).^2 + Z(7,:).^2).';
    g_thr_lo  = d.Tmin - Tmag;
    g_thr_hi  = Tmag - d.Tmax;
    tt        = tan(d.theta_mx);
    g_gs_pos  = ( Z(1,:).' - tt*Z(2,:).');
    g_gs_neg  = (-Z(1,:).' - tt*Z(2,:).');
    c_ineq    = [g_thr_lo; g_thr_hi; g_gs_pos; g_gs_neg];
end

function [t, X] = fwd_integrate(sol, d, mode)
    % Replay the control through ode45, sample at grid nodes.
    %   INPUT
    %     sol  - struct: non-dim solution
    %     d    - struct: non-dim data
    %     mode - 'zoh' (piecewise-constant) or 'pwl' (piecewise-linear)
    %   OUTPUT
    %     t    - node times (Nx1)
    %     X    - replayed state (Nx5)
    arguments
        sol  (1,1) struct
        d    (1,1) struct
        mode {mustBeTextScalar}
    end
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
        rhs_t = @(tt, x) ode_descent(x, u_fcn(tt), d.Vc);
        [~, Y] = ode45(rhs_t, [t_k, t_kp], X(k,:).', opts);
        X(k+1,:) = Y(end,:);
    end
    t = sol.t;
end

function e = node_err(sol, X)
    % Per-node pos+vel error norm, NLP vs ode45 replay (mass excluded).
    %   INPUT
    %     sol - struct: non-dim solution
    %     X   - replayed state (Nx5)
    %   OUTPUT
    %     e   - error norm per node (Nx1)
    e = vecnorm( ...
        [sol.x sol.y sol.vx sol.vy] - X(:,1:4), 2, 2);
end

function plot_compare3(s_t, s_z, s_s, err_t, err_z, err_s, conv_hist, d)
    % Compare trapezoidal / ZOH-RK4 / ZOH-SCvx solutions.
    %   INPUT
    %     s_t,s_z,s_s       - SI solutions (trap, ZOH-RK4, ZOH-SCvx)
    %     err_t,err_z,err_s - per-node fidelity errors
    %     conv_hist         - struct: SCvx trace
    %     d                 - struct: SI problem data
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
    valid = ~isnan(conv_hist.delta);
    iters = find(valid);
    figure('Name','SCvx convergence','Position',[100 100 800 500]);
    tiledlayout(2, 1);
    nexttile;
    yyaxis left;
    semilogy(iters, max(conv_hist.delta(valid), eps), '-o', 'LineWidth', 1.4, 'DisplayName','||\Delta x||');
    hold on;
    semilogy(iters, max(conv_hist.rho(valid),   eps), '-s', 'LineWidth', 1.4, 'DisplayName','\rho (trust scale)');
    ylabel('log scale');
    yyaxis right;
    plot(iters, conv_hist.m_f(valid), '-d', 'LineWidth', 1.4, 'Color',[0.4 0.4 0.4], 'DisplayName','m_f (nondim)');
    ylabel('m_f / m_0');
    grid on;
    legend('Location','best');
    title('SCvx convergence: state delta, trust scale, final mass');

    nexttile;
    bar_h = bar(iters, conv_hist.eta(valid), 0.7);
    ylim([-0.5 2]);
    yline(0.25, 'k--', '\eta_l = 0.25');
    yline(0.7,  'k--', '\eta_h = 0.7');
    if isfield(conv_hist, 'acc')
        for k = 1:length(iters)
            if ~conv_hist.acc(iters(k))
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
% One SCvx inner SOCP via YALMIP + ECOS.
%   INPUT
%     tf         - flight time (nondim)
%     N          - node count
%     d          - struct: non-dim data
%     Abar,Bbar,cbar - discrete LTV matrices
%     ref_sol    - struct: trust-region centre, optional
%     trust      - struct: trust radii, optional
%   OUTPUT
%     sol        - struct: non-dim solution
% LTV dynamics are linear equalities; thrust bound is a SOC; glide-slope
% and trust region are linear.
    arguments
        tf   (1,1) double {mustBePositive, mustBeFinite}
        N    (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(N, 2)}
        d    (1,1) struct
        Abar (5,5,:) double
        Bbar (5,2,:) double
        cbar (5,:)   double
        ref_sol = []
        trust   = []
    end
    X = sdpvar(5, N,   'full');   % state  [x; y; vx; vy; m] at each node
    U = sdpvar(2, N-1, 'full');   % control [Tx; Ty] per ZOH interval

    tt = tan(d.theta_mx);

    cstr = (X(:,1) == [d.x0; d.y0; d.vx0; d.vy0; d.m0]);   % I.C.
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

    if ~isempty(ref_sol) && ~isempty(trust)
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

function [sol, conv_hist] = solve_scvx_yalmip(tf, N, d, max_iter, tol, init_ref, base_trust)
% SCvx outer loop, SOCP inner subproblem (YALMIP/ECOS).
%   INPUT
%     tf         - flight time (nondim)
%     N          - node count
%     d          - struct: non-dim data
%     max_iter   - iteration cap
%     tol        - convergence tol on ||delta x||
%     init_ref   - struct: initial reference, optional
%     base_trust - struct: base trust radii
%   OUTPUT
%     sol        - struct: best non-dim solution (+ iter count)
%     conv_hist  - struct: per-iteration trace
% Same adaptive trust-region logic as solve_scvx; inner solve via
% solve_ltv_nlp_yalmip instead of fmincon.
    arguments
        tf       (1,1) double {mustBePositive, mustBeFinite}
        N        (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(N, 2)}
        d        (1,1) struct
        max_iter (1,1) double {mustBeInteger, mustBePositive}
        tol      (1,1) double {mustBePositive}
        init_ref         = []
        base_trust (1,1) struct = struct('pos', 0.17, 'vel', 0.6, 'mass', 0.1, 'thrust', 1.0)
    end
    idx = @(i) (i-1)*7 + (1:7);
    if isempty(init_ref)
        ref = unpack(init_guess(N, d, idx, true), tf, N, d);
    else
        ref = init_ref;
    end

    rho = 1.0;   rho_min = 1e-3;   rho_max = 1.0;
    eta_l = 0.25;   eta_h = 0.7;

    conv_hist.m_f   = nan(max_iter, 1);   conv_hist.delta = nan(max_iter, 1);
    conv_hist.rho   = nan(max_iter, 1);   conv_hist.eta   = nan(max_iter, 1);
    conv_hist.acc   = false(max_iter, 1);

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
        conv_hist.m_f(iter) = sol_cand.m_f;   conv_hist.delta(iter) = delta_x;
        conv_hist.rho(iter) = rho;             conv_hist.eta(iter)   = eta;
        conv_hist.acc(iter) = accepted;
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

function sol = solve_gfold_socp(tf, N, d, Abar, Bbar, cbar, z_ref, ref_sol, trust)
% One SCvx inner SOCP for the GFOLD log-mass transcription (YALMIP + ECOS).
%   State    XI = [x; y; vx; vy; z],  z = ln(m)
%   Control  W  = [ux; uy; sigma],    u = T/m,  sigma >= ||u||
%   INPUT
%     tf, N, d        - flight time (nondim), node count, non-dim data
%     Abar,Bbar,cbar  - constant LTI ZOH matrices (5x5, 5x3, 5x1)
%     z_ref           - reference log-mass (Nx1) for the thrust-bound linearisation
%     ref_sol, trust  - trust-region centre and radii, optional
%   OUTPUT
%     sol             - struct: non-dim solution (original variables m, T)
% LTI dynamics and the lossless cone ||u||<=sigma are exact; only the upper
% thrust bound sigma <= Tmax*e^{-z} is linearised about z_ref.
    arguments
        tf   (1,1) double {mustBePositive, mustBeFinite}
        N    (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(N, 2)}
        d    (1,1) struct
        Abar (5,5) double
        Bbar (5,3) double
        cbar (5,1) double
        z_ref (:,1) double
        ref_sol = []
        trust   = []
    end
    XI = sdpvar(5, N,   'full');    % [x; y; vx; vy; z]
    W  = sdpvar(3, N-1, 'full');    % [ux; uy; sigma] per ZOH interval
    tt = tan(d.theta_mx);
    z0 = log(d.m0);                 % = 0 since m0_nd = 1

    cstr = (XI(:,1) == [d.x0; d.y0; d.vx0; d.vy0; z0]);   % I.C.
    cstr = [cstr, XI(1:4, N) == 0];                        % terminal BCs (z_N free)

    for k = 1:N-1
        cstr = [cstr, XI(:,k+1) == Abar*XI(:,k) + Bbar*W(:,k) + cbar];   % LTI dynamics
        cstr = [cstr, norm(W(1:2,k)) <= W(3,k)];                          % ||u|| <= sigma (SOC)
        ezr  = exp(-z_ref(k));                                            % linearised thrust upper bound:
        cstr = [cstr, W(3,k) <= d.Tmax*ezr*(1 - (XI(5,k) - z_ref(k)))];  %   sigma <= Tmax*e^{-z_ref}(1-(z-z_ref))
    end

    for k = 1:N     % glide-slope cone (linear)
        cstr = [cstr,  XI(1,k) <= tt*XI(2,k),  -XI(1,k) <= tt*XI(2,k)];
    end
    cstr = [cstr, XI(2,:) >= 0];                          % altitude >= 0
    cstr = [cstr, XI(5,:) >= log(1e-3), XI(5,:) <= 0];    % mass bounds in z = ln(m)

    if ~isempty(ref_sol) && ~isempty(trust)
        for k = 1:N
            cstr = [cstr, XI(1,k) >= ref_sol.x(k)  - trust.pos, XI(1,k) <= ref_sol.x(k)  + trust.pos];
            cstr = [cstr, XI(2,k) >= ref_sol.y(k)  - trust.pos, XI(2,k) <= ref_sol.y(k)  + trust.pos];
            cstr = [cstr, XI(3,k) >= ref_sol.vx(k) - trust.vel, XI(3,k) <= ref_sol.vx(k) + trust.vel];
            cstr = [cstr, XI(4,k) >= ref_sol.vy(k) - trust.vel, XI(4,k) <= ref_sol.vy(k) + trust.vel];
            cstr = [cstr, XI(5,k) >= ref_sol.z(k)  - trust.lz,  XI(5,k) <= ref_sol.z(k)  + trust.lz];
        end
        for k = 1:N-1
            cstr = [cstr, W(1,k) >= ref_sol.ux(k) - trust.u, W(1,k) <= ref_sol.ux(k) + trust.u];
            cstr = [cstr, W(2,k) >= ref_sol.uy(k) - trust.u, W(2,k) <= ref_sol.uy(k) + trust.u];
            cstr = [cstr, W(3,k) >= max(0, ref_sol.sig(k) - trust.sig), W(3,k) <= ref_sol.sig(k) + trust.sig];
        end
    end

    res = optimize(cstr, -XI(5,N), sdpsettings('solver','ecos','verbose',0));
    if res.problem ~= 0 && res.problem ~= 1
        warning('GFOLD YALMIP/ECOS: flag %d (%s)', res.problem, res.info);
    end

    Xv = value(XI);   Wv = value(W);
    sol.t  = linspace(0, tf, N).';
    sol.x  = Xv(1,:).';  sol.y = Xv(2,:).';  sol.vx = Xv(3,:).'; sol.vy = Xv(4,:).';
    sol.z  = Xv(5,:).';  sol.m = exp(sol.z);
    sol.ux = [Wv(1,:).'; 0];  sol.uy = [Wv(2,:).'; 0];  sol.sig = [Wv(3,:).'; 0];
    sol.Tx = sol.m .* sol.ux; sol.Ty = sol.m .* sol.uy;          % T = m*u (last node padded 0)
    sol.Tmag = sqrt(sol.Tx.^2 + sol.Ty.^2);
    sol.tf = tf;   sol.m_f = sol.m(end);   sol.m0 = d.m0;
end

function [sol, conv_hist] = solve_gfold_scvx(tf, N, d, max_iter, tol)
% SCvx outer loop for the GFOLD log-mass transcription.
%   INPUT
%     tf, N, d  - flight time (nondim), node count, non-dim data
%     max_iter  - iteration cap
%     tol       - convergence tol on ||delta xi||
%   OUTPUT
%     sol       - struct: best non-dim solution (+ iter count)
%     conv_hist - struct: per-iteration trace
% Dynamics are exact LTI (lti_zoh, computed once), so the loop linearises
% ONLY the upper thrust bound about the current z and self-starts from the
% analytic max-thrust mass profile (no trapezoidal warm start needed). The
% adaptive trust region / ratio test mirror solve_scvx_yalmip.
    arguments
        tf       (1,1) double {mustBePositive, mustBeFinite}
        N        (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(N, 2)}
        d        (1,1) struct
        max_iter (1,1) double {mustBeInteger, mustBePositive}
        tol      (1,1) double {mustBePositive}
    end
    [Abar, Bbar, cbar] = lti_zoh(tf/(N-1), d.Vc);      % exact LTI ZOH, once
    t_grid = linspace(0, tf, N).';

    % Self-starting reference: analytic max-thrust mass-depletion profile for
    % z, linear interpolation for position/velocity, hover-ish acceleration.
    m_apri = max(d.m0 - d.Vc*d.Tmax*t_grid, 1e-2);
    ref.z  = log(m_apri);
    ref.x  = linspace(d.x0, 0, N).';   ref.y  = linspace(d.y0, 0, N).';
    ref.vx = linspace(d.vx0, 0, N).';  ref.vy = linspace(d.vy0, 0, N).';
    ref.m  = exp(ref.z);
    ref.ux = zeros(N,1);  ref.uy = ones(N,1);  ref.sig = ones(N,1);
    ref.m_f = exp(ref.z(N));

    base = struct('pos', 0.5, 'vel', 1.0, 'lz', 0.4, 'u', 4.0, 'sig', 4.0);
    rho = 1.0;   rho_min = 1e-3;   rho_max = 1.0;   eta_l = 0.25;   eta_h = 0.7;

    conv_hist.m_f   = nan(max_iter, 1);   conv_hist.delta = nan(max_iter, 1);
    conv_hist.rho   = nan(max_iter, 1);   conv_hist.eta   = nan(max_iter, 1);
    conv_hist.acc   = false(max_iter, 1);

    sol_best = ref;   converged = false;
    for iter = 1:max_iter
        if iter == 1
            % First solve free of the trust region: dynamics are exact, so let
            % the SOCP find a dynamically feasible trajectory before refining.
            cand = solve_gfold_socp(tf, N, d, Abar, Bbar, cbar, ref.z, [], []);
        else
            scaled = struct('pos', rho*base.pos, 'vel', rho*base.vel, 'lz', rho*base.lz, ...
                            'u', rho*base.u, 'sig', rho*base.sig);
            cand = solve_gfold_socp(tf, N, d, Abar, Bbar, cbar, ref.z, ref, scaled);
        end

        J_pred = cand.m_f - ref.m_f;
        [~, X_act] = fwd_integrate_uacc(cand, d);        % nonlinear u-ZOH replay
        J_act  = X_act(end, 5) - ref.m_f;
        eta    = ternary(abs(J_pred) < 1e-10, 1, J_act / J_pred);

        delta_x = norm([cand.x - ref.x; cand.y - ref.y; cand.vx - ref.vx; ...
                        cand.vy - ref.vy; cand.z - ref.z]);

        accepted = (eta >= eta_l);
        conv_hist.m_f(iter) = cand.m_f;   conv_hist.delta(iter) = delta_x;
        conv_hist.rho(iter) = rho;         conv_hist.eta(iter)   = eta;
        conv_hist.acc(iter) = accepted;
        fprintf('  SCvx-GFOLD iter %2d:  rho=%.3f  eta=%+7.3f  delta_x=%.3e  m_f=%.4f  %s\n', ...
            iter, rho, eta, delta_x, cand.m_f, ternary(accepted,'ACCEPTED','rejected'));

        if accepted
            sol_best = cand;   ref = cand;   ref.m_f = cand.m_f;
            if eta > eta_h,  rho = min(rho_max, 2*rho);  end
            if delta_x < tol
                fprintf('    SCvx-GFOLD converged (delta_x < tol).\n');
                converged = true;   break;
            end
        else
            rho = 0.5 * rho;
            if rho < rho_min
                fprintf('    SCvx-GFOLD stopping: trust region collapsed.\n');   break;
            end
        end
    end
    if ~converged && iter == max_iter
        fprintf('    SCvx-GFOLD hit iteration cap (no further improvement).\n');
    end
    sol = sol_best;   sol.iter = iter;
end

function [t, X] = fwd_integrate_uacc(sol, d)
    % Replay holding the acceleration u = T/m piecewise constant (T = m(t)*u
    % floats) -- the ZOH convention native to GFOLD. Mirrors fwd_integrate,
    % but the control is the acceleration and the RHS is ode_descent_uacc.
    arguments
        sol (1,1) struct
        d   (1,1) struct
    end
    N = numel(sol.t);   X = zeros(N, 5);
    X(1,:) = [d.x0, d.y0, d.vx0, d.vy0, d.m0];
    opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
    for k = 1:N-1
        uacc   = [sol.ux(k); sol.uy(k)];
        rhs_t  = @(tt, x) ode_descent_uacc(x, uacc, d.Vc);
        [~, Y] = ode45(rhs_t, [sol.t(k), sol.t(k+1)], X(k,:).', opts);
        X(k+1,:) = Y(end,:);
    end
    t = sol.t;
end

function plot_compare5(s_t, s_z, s_s, s_y, s_g, err_t, err_z, err_s, err_y, err_g, hist_s, hist_y, hist_g, d)
    % Same as plot_compare3 plus the YALMIP/ECOS and GFOLD log-mass variants.
    %   INPUT
    %     s_t,s_z,s_s,s_y,s_g           - SI solutions (trap, ZOH-RK4, SCvx,
    %                                     SCvx-YALMIP, GFOLD log-mass)
    %     err_t,err_z,err_s,err_y,err_g - per-node fidelity errors
    %     hist_s,hist_y,hist_g          - SCvx traces (fmincon, YALMIP, GFOLD)
    %     d                             - struct: SI problem data
    cT = [0.0  0.4  0.8];
    cZ = [0.85 0.33 0.1];
    cS = [0.47 0.67 0.19];
    cY = [0.49 0.18 0.56];   % purple — YALMIP variant
    cG = [0.0  0.6  0.55];   % teal — GFOLD log-mass variant

    % --- Trajectory ---
    figure('Name','Trajectory comparison','Position',[100 100 600 500]);
    hold on; grid on; axis equal;
    yy = linspace(0, max([s_t.y; s_z.y; s_s.y; s_y.y; s_g.y])*1.05, 50);
    xx = tan(d.theta_mx) * yy;
    plot( xx, yy, '--', 'Color',[0.5 0.5 0.5], 'HandleVisibility','off');
    plot(-xx, yy, '--', 'Color',[0.5 0.5 0.5], 'HandleVisibility','off');
    plot(s_t.x, s_t.y, '-',   'Color', cT, 'LineWidth', 1.6, 'DisplayName','Trapezoidal');
    plot(s_z.x, s_z.y, '--',  'Color', cZ, 'LineWidth', 1.6, 'DisplayName','ZOH (RK4)');
    plot(s_s.x, s_s.y, ':',   'Color', cS, 'LineWidth', 2.0, 'DisplayName','ZOH (LTV+SCvx)');
    plot(s_y.x, s_y.y, '-.',  'Color', cY, 'LineWidth', 2.0, 'DisplayName','ZOH (LTV+SCvx YALMIP)');
    plot(s_g.x, s_g.y, '-',   'Color', cG, 'LineWidth', 1.4, 'DisplayName','GFOLD (log-mass)');
    plot(0, 0, 'k^', 'MarkerSize', 8, 'MarkerFaceColor','k', 'HandleVisibility','off');
    xlabel('x  [m]'); ylabel('y  [m]');
    title('Descent trajectory: five transcriptions vs glide-slope corridor');
    legend('Location','best');

    % --- Thrust magnitude ---
    figure('Name','Thrust comparison','Position',[100 100 600 400]);
    hold on; grid on;
    plot(  s_t.t, s_t.Tmag/1e3, '-',  'Color', cT, 'LineWidth', 1.6, 'DisplayName','Trapezoidal (PWL)');
    stairs(s_z.t, s_z.Tmag/1e3, '--', 'Color', cZ, 'LineWidth', 1.6, 'DisplayName','ZOH RK4 (PWC)');
    stairs(s_s.t, s_s.Tmag/1e3, ':',  'Color', cS, 'LineWidth', 2.0, 'DisplayName','ZOH LTV+SCvx (PWC)');
    stairs(s_y.t, s_y.Tmag/1e3, '-.', 'Color', cY, 'LineWidth', 2.0, 'DisplayName','ZOH LTV+SCvx YALMIP (PWC)');
    stairs(s_g.t, s_g.Tmag/1e3, '-',  'Color', cG, 'LineWidth', 1.4, 'DisplayName','GFOLD log-mass (PWC)');
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
    plot(s_g.t, s_g.m, '-',  'Color', cG, 'LineWidth', 1.4, 'DisplayName','GFOLD (log-mass)');
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
    semilogy(s_g.t, max(err_g, eps), '-',  'Color', cG, 'LineWidth', 1.4, ...
        'DisplayName','GFOLD log-mass (u-ZOH fwd-int)');
    xlabel('t  [s]'); ylabel('||x_{NLP} - x_{ode45}|| (nondim)');
    title('Transcription fidelity at grid nodes (non-dim state norm)');
    legend('Location','best');

    % --- SCvx convergence: fmincon vs YALMIP vs GFOLD log-mass ---
    for ii = 1:3
        switch ii
            case 1,  conv_hist = hist_s;  lbl = 'fmincon/SQP';     col = cS;
            case 2,  conv_hist = hist_y;  lbl = 'YALMIP/ECOS';     col = cY;
            case 3,  conv_hist = hist_g;  lbl = 'GFOLD log-mass';  col = cG;
        end
        valid = ~isnan(conv_hist.delta);
        iters = find(valid);
        figure('Name', sprintf('SCvx convergence — %s', lbl), 'Position',[100+400*(ii-1) 600 800 500]);
        tiledlayout(2, 1);
        nexttile;
        yyaxis left;
        semilogy(iters, max(conv_hist.delta(valid), eps), '-o', 'LineWidth', 1.4, ...
            'Color', col, 'DisplayName','||\Delta x||');
        hold on;
        semilogy(iters, max(conv_hist.rho(valid), eps), '-s', 'LineWidth', 1.4, ...
            'Color', col*0.6, 'DisplayName','\rho (trust scale)');
        ylabel('log scale');
        yyaxis right;
        plot(iters, conv_hist.m_f(valid), '-d', 'LineWidth', 1.4, 'Color',[0.4 0.4 0.4], ...
            'DisplayName','m_f (nondim)');
        ylabel('m_f / m_0');
        grid on; legend('Location','best');
        title(sprintf('SCvx convergence [%s]: state delta, trust scale, final mass', lbl));
        nexttile;
        bar_h = bar(iters, conv_hist.eta(valid), 0.7);
        bar_h.FaceColor = col;
        ylim([-0.5 2]);
        yline(0.25, 'k--', '\eta_l = 0.25');
        yline(0.7,  'k--', '\eta_h = 0.7');
        if isfield(conv_hist, 'acc')
            for k = 1:length(iters)
                if ~conv_hist.acc(iters(k))
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
