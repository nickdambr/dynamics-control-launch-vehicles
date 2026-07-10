%% PROTOTYPE — GFOLD powered descent via log-mass change of variables + SCvx
%
% NOTE: this exploratory prototype has been promoted to an official Task-2
% "variant (d)" inside main_task2.m (solve_gfold_scvx / solve_gfold_socp), and
% the numeric kernels now live in standalone files lti_zoh.m and
% ode_descent_uacc.m (tested by tests/gfoldLogMassTest.m). The script is kept
% as a minimal self-contained GFOLD demo that runs in a few seconds without the
% other four transcriptions; it reuses those standalone kernels. Safe to delete
% if the integrated version is enough.
%
% Change of variables (GFOLD, Acikmese/Blackmore):
%   state    xi = [x; y; vx; vy; z],  z = ln(m)
%   control  w  = [ux; uy; sigma],    u = T/m (acceleration), sigma >= ||u||
%
% Under this substitution the non-dim dynamics become LINEAR TIME-INVARIANT:
%   xdot = vx ;  ydot = vy ;  vxdot = ux ;  vydot = uy - 1 ;  zdot = -Vc*sigma
% so the ZOH discretisation is exact and computed ONCE (van Loan / expm),
% removing the epsilon-regularised singular mass row of variants (b)/(c).
%
% The ONLY residual non-convexity is the upper thrust bound
%   ||u_k|| <= Tmax * exp(-z_k)        (region below a convex curve -> non-convex)
% which is linearised about a reference z_ref and iterated in an SCvx outer
% loop with an adaptive box trust region and an ode45 ratio test.
%
% Reuses ode_descent.m only as a "ground-truth" cross-check.
%
% Niccolo D'Ambrosio — HM2 prototype.

clear; clc;

%% Problem data (Table 1) and non-dimensionalisation
d_si.x0   = 1000;   d_si.y0   = 3000;
d_si.vx0  = 300;    d_si.vy0  = -200;
d_si.m0   = 2000;
d_si.g    = 9.81;
d_si.Isp  = 225;    d_si.g0   = 9.80665;
d_si.c    = d_si.Isp * d_si.g0;
d_si.Tmin = 0;      d_si.Tmax = 70000;
d_si.theta_mx = deg2rad(60);
tf_si = 38;         N = 50;

[ref, d] = nondim(d_si);            % ref scales + non-dim data struct
tf  = tf_si / ref.t;                % non-dim flight time
dt  = tf / (N - 1);

fprintf('=== GFOLD log-mass prototype ===\n');
fprintf('Non-dim: L=%.1f m  V=%.2f m/s  t=%.3f s  Vc=%.4f\n', ref.L, ref.V, ref.t, d.Vc);
fprintf('Tmax_nd=%.4f  tf_nd=%.4f  dt_nd=%.4f  N=%d\n\n', d.Tmax, tf, dt, N);

%% Solver availability
if ~(exist('sdpvar','file') && exist('ecos','file'))
    error('YALMIP and/or ECOS not found on the path — cannot run the SOCP prototype.');
end

%% SCvx solve (the exact LTI ZOH matrices are built once inside the solver)
max_iter = 20;   tol = 1e-3;
sol = solve_gfold_scvx(tf, N, d, max_iter, tol);

%% Validation in ORIGINAL variables
sol_dim = dim_sol(sol, ref);
[~, X_rep] = fwd_integrate_uacc(sol, d);            % u-ZOH nonlinear replay (non-dim)
pos_err = norm(X_rep(end,1:2)) * ref.L;             % touchdown position error [m]
vel_err = norm(X_rep(end,3:4)) * ref.V;             % touchdown velocity error [m/s]
dmf     = (sol.m_f - X_rep(end,5)) * ref.m;         % model vs replay final-mass drift [kg]
node_e  = vecnorm([sol.x sol.y sol.vx sol.vy] - X_rep(:,1:4), 2, 2);

% Constraint checks (original variables, SI)
Tmag_si   = sol_dim.Tmag;                            % thrust magnitude per node [N]
gs_margin = tan(d.theta_mx)*sol.y - abs(sol.x);      % >=0 means glide-slope satisfied
mass_mono = all(diff(sol.m) <= 1e-9);                % mass non-increasing
cone_gap  = max(abs(sol.sig(1:N-1) - vecnorm([sol.ux(1:N-1) sol.uy(1:N-1)],2,2)));

fprintf('\n=== RESULTS ===\n');
fprintf('m_f (model)        : %.2f kg   (fuel %.2f kg)\n', sol_dim.m_f, sol_dim.fuel);
fprintf('m_f (ode45 replay) : %.2f kg   (drift %.3f kg)\n', X_rep(end,5)*ref.m, dmf);
fprintf('Reference (Task 2) : trap 1403.2  ZOH-RK4 1403.4  SCvx 1399.5-1400.8 kg\n');
fprintf('Touchdown replay   : pos %.4f m   vel %.5f m/s\n', pos_err, vel_err);
fprintf('Max node error     : %.3e (non-dim pos+vel norm)\n', max(node_e));
fprintf('Max thrust         : %.1f N  (limit %.0f N) -> %s\n', max(Tmag_si), d_si.Tmax, ...
        ternary(max(Tmag_si) <= d_si.Tmax + 1, 'OK', 'VIOLATED'));
fprintf('Glide-slope        : min margin %.4e (>=0 OK) -> %s\n', min(gs_margin), ...
        ternary(min(gs_margin) >= -1e-6, 'OK', 'VIOLATED'));
fprintf('Altitude y>=0      : min y %.4e -> %s\n', min(sol.y), ...
        ternary(min(sol.y) >= -1e-6, 'OK', 'VIOLATED'));
fprintf('Mass monotone      : %s\n', ternary(mass_mono, 'OK', 'VIOLATED'));
fprintf('Lossless cone gap  : %.3e (||u||=sigma at optimum)\n', cone_gap);

%% Quick figure
figure('Name','GFOLD log-mass prototype','Color','w');
subplot(2,2,1); plot(sol_dim.x, sol_dim.y, '-o','MarkerSize',3); grid on;
xlabel('x [m]'); ylabel('y [m]'); title('Trajectory'); axis equal;
subplot(2,2,2); stairs(sol_dim.t, sol_dim.Tmag/1e3, 'LineWidth',1.4); grid on;
xlabel('t [s]'); ylabel('||T|| [kN]'); title('Thrust (PWC)'); yline(d_si.Tmax/1e3,'r--');
subplot(2,2,3); plot(sol_dim.t, sol_dim.m, 'LineWidth',1.4); grid on;
xlabel('t [s]'); ylabel('m [kg]'); title('Mass');
subplot(2,2,4); semilogy(sol_dim.t, max(node_e,eps), 'LineWidth',1.4); grid on;
xlabel('t [s]'); ylabel('node error'); title('Transcription fidelity (u-ZOH replay)');

% ===================================================================== %
%                           LOCAL FUNCTIONS                              %
% ===================================================================== %

function [ref, dnd] = nondim(d)
    % Reference scales + non-dim data (mirrors main_task2.m:nondim).
    ref.L = d.y0;            ref.g = d.g;
    ref.t = sqrt(ref.L/ref.g);   ref.V = sqrt(ref.g*ref.L);
    ref.m = d.m0;            ref.T = ref.m*ref.g;
    dnd.x0  = d.x0/ref.L;    dnd.y0  = d.y0/ref.L;
    dnd.vx0 = d.vx0/ref.V;   dnd.vy0 = d.vy0/ref.V;
    dnd.m0  = d.m0/ref.m;    dnd.Tmin = d.Tmin/ref.T;
    dnd.Tmax = d.Tmax/ref.T; dnd.Vc = ref.V/d.c;
    dnd.theta_mx = d.theta_mx;
end

function sol = dim_sol(s, ref)
    % Scale a non-dim GFOLD solution back to SI (mirrors main_task2.m:dim_sol).
    sol.t  = s.t*ref.t;     sol.x = s.x*ref.L;   sol.y = s.y*ref.L;
    sol.vx = s.vx*ref.V;    sol.vy = s.vy*ref.V; sol.m = s.m*ref.m;
    sol.Tx = s.Tx*ref.T;    sol.Ty = s.Ty*ref.T;
    sol.Tmag = sqrt(sol.Tx.^2 + sol.Ty.^2);
    sol.tf = s.tf*ref.t;    sol.m_f = s.m_f*ref.m;
    sol.fuel = (s.m0 - s.m_f)*ref.m;
end

function sol = solve_gfold_socp(tf, N, d, Abar, Bbar, cbar, z_ref, ref_sol, trust)
    % One SCvx inner SOCP (YALMIP+ECOS) in log-mass GFOLD variables.
    XI = sdpvar(5, N,   'full');    % [x; y; vx; vy; z]
    W  = sdpvar(3, N-1, 'full');    % [ux; uy; sigma] per ZOH interval
    tt = tan(d.theta_mx);
    z0 = log(d.m0);                 % = 0 since m0_nd = 1

    cstr = (XI(:,1) == [d.x0; d.y0; d.vx0; d.vy0; z0]);   % I.C.
    cstr = [cstr, XI(1:4,N) == 0];                         % terminal BCs (z_N free)

    for k = 1:N-1
        cstr = [cstr, XI(:,k+1) == Abar*XI(:,k) + Bbar*W(:,k) + cbar];   % LTI dynamics
        cstr = [cstr, norm(W(1:2,k)) <= W(3,k)];                          % ||u|| <= sigma (SOC)
        % linearised upper thrust bound: sigma <= Tmax*e^{-z_ref}(1-(z-z_ref))
        ezr  = exp(-z_ref(k));
        cstr = [cstr, W(3,k) <= d.Tmax*ezr*(1 - (XI(5,k) - z_ref(k)))];
    end

    for k = 1:N                                                          % glide-slope (linear)
        cstr = [cstr,  XI(1,k) <= tt*XI(2,k),  -XI(1,k) <= tt*XI(2,k)];
    end
    cstr = [cstr, XI(2,:) >= 0];                                         % altitude
    cstr = [cstr, XI(5,:) >= log(1e-3), XI(5,:) <= 0];                   % mass (z) bounds

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
            cstr = [cstr, W(3,k) >= max(0,ref_sol.sig(k) - trust.sig), W(3,k) <= ref_sol.sig(k) + trust.sig];
        end
    end

    res = optimize(cstr, -XI(5,N), sdpsettings('solver','ecos','verbose',0));
    if res.problem ~= 0
        warning('YALMIP/ECOS: flag %d (%s)', res.problem, res.info);
    end

    Xv = value(XI);   Wv = value(W);
    sol.t  = linspace(0, tf, N).';
    sol.x  = Xv(1,:).';  sol.y = Xv(2,:).';  sol.vx = Xv(3,:).'; sol.vy = Xv(4,:).';
    sol.z  = Xv(5,:).';  sol.m = exp(sol.z);
    sol.ux = [Wv(1,:).'; 0];  sol.uy = [Wv(2,:).'; 0];  sol.sig = [Wv(3,:).'; 0];
    sol.Tx = sol.m.*sol.ux;   sol.Ty = sol.m.*sol.uy;    % T = m*u (last node padded 0)
    sol.Tmag = sqrt(sol.Tx.^2 + sol.Ty.^2);
    sol.tf = tf;   sol.m_f = sol.m(end);   sol.m0 = d.m0;
end

function [sol, hist] = solve_gfold_scvx(tf, N, d, max_iter, tol)
    % SCvx outer loop around the GFOLD SOCP (mirrors solve_scvx_yalmip).
    [Abar, Bbar, cbar] = lti_zoh(tf/(N-1), d.Vc);   % exact LTI ZOH, once
    t_grid = linspace(0, tf, N).';

    % --- warm-start reference ---
    % z_ref0: a-priori max-thrust mass-depletion profile (classic GFOLD ref).
    m_apri = max(d.m0 - d.Vc*d.Tmax*t_grid, 1e-2);
    ref.z  = log(m_apri);
    ref.x  = linspace(d.x0, 0, N).';   ref.y  = linspace(d.y0, 0, N).';
    ref.vx = linspace(d.vx0, 0, N).';  ref.vy = linspace(d.vy0, 0, N).';
    ref.m  = exp(ref.z);
    ref.ux = zeros(N,1);  ref.uy = ones(N,1);  ref.sig = ones(N,1);   % hover-ish accel
    ref.m_f = exp(ref.z(N));

    base = struct('pos',0.5, 'vel',1.0, 'lz',0.4, 'u',4.0, 'sig',4.0);
    rho = 1.0;  rho_min = 1e-3;  rho_max = 1.0;  eta_l = 0.25;  eta_h = 0.7;

    hist.m_f = nan(max_iter,1); hist.delta = nan(max_iter,1);
    hist.rho = nan(max_iter,1); hist.eta = nan(max_iter,1); hist.acc = false(max_iter,1);

    sol_best = ref;  converged = false;
    for it = 1:max_iter
        if it == 1
            % First solve free of the trust region: the dynamics are exact, so
            % let the SOCP find a dynamically feasible trajectory before the
            % trust-region refinement kicks in (avoids infeasibility against a
            % crude, non-dynamic warm-start).
            cand = solve_gfold_socp(tf, N, d, Abar, Bbar, cbar, ref.z, [], []);
        else
            tr = struct('pos',rho*base.pos,'vel',rho*base.vel,'lz',rho*base.lz, ...
                        'u',rho*base.u,'sig',rho*base.sig);
            cand = solve_gfold_socp(tf, N, d, Abar, Bbar, cbar, ref.z, ref, tr);
        end

        J_pred = cand.m_f - ref.m_f;
        [~, X_act] = fwd_integrate_uacc(cand, d);          % nonlinear u-ZOH replay
        J_act  = X_act(end,5) - ref.m_f;
        eta    = ternary(abs(J_pred) < 1e-10, 1, J_act/J_pred);

        delta = norm([cand.x-ref.x; cand.y-ref.y; cand.vx-ref.vx; cand.vy-ref.vy; cand.z-ref.z]);

        accepted = (eta >= eta_l);
        hist.m_f(it)=cand.m_f; hist.delta(it)=delta; hist.rho(it)=rho; hist.eta(it)=eta; hist.acc(it)=accepted;
        fprintf('  SCvx-GFOLD it %2d:  rho=%.3f  eta=%+8.3f  delta=%.3e  m_f=%.4f  %s\n', ...
            it, rho, eta, delta, cand.m_f, ternary(accepted,'ACCEPTED','rejected'));

        if accepted
            sol_best = cand;  ref = cand;  ref.m_f = cand.m_f;
            if eta > eta_h, rho = min(rho_max, 2*rho); end
            if delta < tol
                fprintf('    converged (delta < tol).\n'); converged = true; break;
            end
        else
            rho = 0.5*rho;
            if rho < rho_min, fprintf('    trust region collapsed.\n'); break; end
        end
    end
    if ~converged && it == max_iter, fprintf('    hit iteration cap.\n'); end
    sol = sol_best;  sol.iter = it;
end

function [t, X] = fwd_integrate_uacc(sol, d)
    % Replay holding the acceleration u constant over each interval (T=m(t)*u
    % floats) — the ZOH convention native to GFOLD. Mirrors fwd_integrate.
    N = numel(sol.t);  X = zeros(N,5);
    X(1,:) = [d.x0, d.y0, d.vx0, d.vy0, d.m0];
    opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
    for k = 1:N-1
        uacc = [sol.ux(k); sol.uy(k)];
        rhs  = @(tt,x) ode_descent_uacc(x, uacc, d.Vc);
        [~, Y] = ode45(rhs, [sol.t(k), sol.t(k+1)], X(k,:).', opts);
        X(k+1,:) = Y(end,:);
    end
    t = sol.t;
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
