%% VALIDATION - Staging corner condition (Appendix C, Step 7)
%  Closes the free staging time t_s with the interior-point corner condition
%  instead of the parametric sweep of main_task4.m. t_s becomes a fifth shooting
%  unknown and the corner condition a fifth residual; the augmented solve must
%  reproduce the swept optimum (ts ~ 0.336).
%
%  Corner condition (report, Task 4):
%     eta*[lam_m(tf) - lam_m(ts)] = c*||lam_v(ts)||*(1/m_plus - 1/m_minus)
%  with  m_minus = 1 - Q*ts,  m_plus = 1 - (1+eta)*Q*ts.
%
%  Reuses ode_burn.m (same RHS as main_task4.m); run from HM1/.

clear; close all; clc;

%% Parameters (identical to main_task4.m)
c   = 0.6;
eta = 0.1;
yf  = 0.04;
Q   = 2;
T   = c * Q;

opts_ode = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
opts_fs  = optimoptions('fsolve', 'Display', 'off', ...
    'MaxIterations', 500, 'MaxFunctionEvaluations', 5000, ...
    'FunctionTolerance', 1e-12, 'StepTolerance', 1e-12);

p.c = c; p.Q = Q; p.T = T; p.yf = yf; p.eta = eta;

%% Step 1: warm-start from a fixed-ts inner BVP (4 unknowns)
ts0 = 0.33;
p_inner = p; p_inner.ts = ts0;
z4_guess = [0.6; 3.8; 14; 0.42];   % [lam_vx0; lam_vy0; lam_y; tf]
[z4, ~, ef1] = fsolve(@(z) shooting_inner(z, p_inner, opts_ode), z4_guess, opts_fs);
if ef1 <= 0
    error('Inner BVP warm-start did not converge.');
end

%% Step 2: augmented solve - 5 unknowns [lam_vx0; lam_vy0; lam_y; tf; ts]
%  Residuals: 3 terminal + H0=0 at t0 (free tf) + corner condition at ts
w_guess = [z4; ts0];
[w, res, ef] = fsolve(@(w) shooting_corner(w, p, opts_ode), w_guess, opts_fs);

%% Extract and report
tf = w(4); ts = w(5);
[zf, ~] = propagate(w, p, opts_ode);
mf = zf(5);
payload = mf - eta*Q*(tf - ts);

fprintf('=== Staging corner-condition solve (5 unknowns, 5 residuals) ===\n');
fprintf('  exitflag        : %d\n', ef);
fprintf('  staging time ts : %.6f   (sweep: 0.336)\n', ts);
fprintf('  burn time    tf : %.6f   (sweep: 0.424)\n', tf);
fprintf('  payload      mu : %.6f   (sweep: 0.068)\n', payload);
fprintf('  residual norm   : %.3e\n', norm(res));

if ef > 0 && abs(ts - 0.336) < 5e-3
    fprintf('\nPASS: corner-condition solve matches the swept optimum.\n');
else
    fprintf('\nFAIL: corner-condition solve did not match the sweep.\n');
end

%% Check: the burnout reference lam_m(tf) is essential.
%  Replacing eta*[lam_m(tf)-lam_m(ts)] by eta*lam_m(ts) alone (the inner-BVP
%  costate, normalized at lam_m0=1) misplaces the optimum.
[w_bad, ~, ef_bad] = fsolve(@(w) shooting_corner_wrong(w, p, opts_ode), [z4; ts0], opts_fs);
fprintf('\n--- Cautionary check: un-referenced (wrong-frame) corner condition ---\n');
if ef_bad > 0
    fprintf('  eta*lam_m(ts) alone  => ts = %.6f  (spurious)\n', w_bad(5));
    fprintf('  correct (burnout-ref) => ts = %.6f\n', ts);
else
    fprintf('  (wrong-frame solve did not converge)\n');
end

%% ===================== LOCAL FUNCTIONS =====================
function [zf, zs_minus] = propagate(w, p, opts_ode)
% Integrate stage 1 [0,ts], jettison structure, stage 2 [ts,tf].
%   INPUT
%     w        - unknowns [lam_vx0; lam_vy0; lam_y; tf; ts]
%     p        - struct: c, Q, T, yf, eta
%     opts_ode - ode45 options
%   OUTPUT
%     zf       - terminal state (6x1)
%     zs_minus - state at ts^- (pre-staging, 6x1)
    pp = p;
    pp.lam_vx0 = w(1); pp.lam_vy0 = w(2); pp.lam_y = w(3);
    tf = w(4); ts = w(5);
    ic = [0;0;0;0;1;1];                       % lam_m0 = 1 (normalization)
    [~, Z1] = ode45(@(t,z) ode_burn(t,z,pp), [0 ts], ic, opts_ode);
    zs_minus = Z1(end,:)';
    z_plus = zs_minus;
    z_plus(5) = z_plus(5) - p.eta*p.Q*ts;     % jettison stage-1 structure
    [~, Z2] = ode45(@(t,z) ode_burn(t,z,pp), [ts tf], z_plus, opts_ode);
    zf = Z2(end,:)';
end

function res = shooting_inner(z, p, opts_ode)
% Fixed-ts inner BVP, same formulation as main_task4.m.
%   INPUT
%     z        - unknowns [lam_vx0; lam_vy0; lam_y; tf]
%     p        - struct: c, Q, T, yf, eta, ts
%     opts_ode - ode45 options
%   OUTPUT
%     res - 4 residuals [y(tf)-yf; vx(tf)-1; vy(tf); H0]
% No arguments block by design: runs inside the fsolve loop.
    lam_vx0 = z(1); lam_vy0 = z(2); lam_y = z(3); tf = z(4); ts = p.ts;
    if tf <= ts || tf > 2 || ts <= 0
        res = 1e6*ones(4,1); return;
    end
    w = [lam_vx0; lam_vy0; lam_y; tf; ts];
    try
        [zf, zs] = propagate(w, p, opts_ode);
    catch
        res = 1e6*ones(4,1); return;
    end
    if zs(5) - p.eta*p.Q*ts <= 0
        res = 1e6*ones(4,1); return;
    end
    Lam0 = sqrt(lam_vx0^2 + lam_vy0^2);
    H0 = -lam_vy0 + p.T*(Lam0 - 1/p.c);       % H=0 at t0 (vx0=vy0=0, m0=1, lam_m0=1)
    res = [zf(2)-p.yf; zf(3)-1; zf(4); H0];
end

function res = shooting_corner(w, p, opts_ode)
% Augmented residual: inner BVP (4) + corner condition (1).
%   INPUT
%     w        - unknowns [lam_vx0; lam_vy0; lam_y; tf; ts]
%     p        - struct: c, Q, T, yf, eta
%     opts_ode - ode45 options
%   OUTPUT
%     res - 5 residuals [y(tf)-yf; vx(tf)-1; vy(tf); H0; corner]
% No arguments block by design: runs inside the fsolve loop.
    lam_vx0 = w(1); lam_vy0 = w(2); lam_y = w(3); tf = w(4); ts = w(5);
    if tf <= ts || tf > 2 || ts <= 0 || ts > 1/(p.Q*(1+p.eta)) - 1e-3
        res = 1e6*ones(5,1); return;
    end
    try
        [zf, zs] = propagate(w, p, opts_ode);
    catch
        res = 1e6*ones(5,1); return;
    end
    m_minus = zs(5);
    m_plus  = m_minus - p.eta*p.Q*ts;
    if m_plus <= 0
        res = 1e6*ones(5,1); return;
    end
    % Terminal conditions + free-tf transversality (H = 0 at t0)
    Lam0 = sqrt(lam_vx0^2 + lam_vy0^2);
    H0   = -lam_vy0 + p.T*(Lam0 - 1/p.c);
    % Corner condition at ts (lam_m continuous, burnout-referenced)
    lam_vy_ts = lam_vy0 - lam_y*ts;
    lam_v_ts  = sqrt(lam_vx0^2 + lam_vy_ts^2);
    corner = p.eta*(zf(6) - zs(6)) - p.c*lam_v_ts*(1/m_plus - 1/m_minus);
    res = [zf(2)-p.yf; zf(3)-1; zf(4); H0; corner];
end

function res = shooting_corner_wrong(w, p, opts_ode)
% Like shooting_corner but with the wrong corner term eta*lam_m(ts) instead of
% eta*[lam_m(tf)-lam_m(ts)]. Shows that dropping the burnout reference misplaces
% the optimum.
%   INPUT
%     w        - unknowns [lam_vx0; lam_vy0; lam_y; tf; ts]
%     p        - struct: c, Q, T, yf, eta
%     opts_ode - ode45 options
%   OUTPUT
%     res - 5 residuals (last one uses the un-referenced corner term)
% No arguments block by design: runs inside the fsolve loop.
    lam_vx0 = w(1); lam_vy0 = w(2); lam_y = w(3); tf = w(4); ts = w(5);
    if tf <= ts || tf > 2 || ts <= 0 || ts > 1/(p.Q*(1+p.eta)) - 1e-3
        res = 1e6*ones(5,1); return;
    end
    try
        [zf, zs] = propagate(w, p, opts_ode);
    catch
        res = 1e6*ones(5,1); return;
    end
    m_minus = zs(5);
    m_plus  = m_minus - p.eta*p.Q*ts;
    if m_plus <= 0
        res = 1e6*ones(5,1); return;
    end
    Lam0 = sqrt(lam_vx0^2 + lam_vy0^2);
    H0   = -lam_vy0 + p.T*(Lam0 - 1/p.c);
    lam_vy_ts = lam_vy0 - lam_y*ts;
    lam_v_ts  = sqrt(lam_vx0^2 + lam_vy_ts^2);
    corner = p.eta*zs(6) - p.c*lam_v_ts*(1/m_plus - 1/m_minus);   % wrong: lam_m(ts) only
    res = [zf(2)-p.yf; zf(3)-1; zf(4); H0; corner];
end
