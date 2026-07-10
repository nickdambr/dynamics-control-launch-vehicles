%% VALIDATION - Rocket sled (shooting benchmark, Appendix A)
%  Min-energy point-to-point transfer with a closed-form solution. Checks the
%  ode45 + fsolve single-shooting machinery (same tolerances as the ascent BVP)
%  against the analytic optimum before using it on the ascent problem.
%
%  Problem:  rdot = v,  vdot = u,   min J = int_0^2 u^2 dt
%            r(0)=0, v(0)=0, r(2)=1/2, v(2)=0,  tf = 2 fixed.
%  Hamiltonian:  H = -u^2 + lam_r*v + lam_v*u
%  Costates:  lam_r_dot = 0,  lam_v_dot = -lam_r ;  u = lam_v/2.
%  Closed form: lam_r0 = lam_v0 = 3/2  =>  u*(t) = (3/4)(1 - t).

clear; close all; clc;

%% Problem data
tf = 2;
rf = 1/2;     % target displacement
vf = 0;       % target velocity

opts_ode = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
opts_fs  = optimoptions('fsolve', 'Display', 'off', ...
    'MaxIterations', 500, 'MaxFunctionEvaluations', 5000, ...
    'FunctionTolerance', 1e-12, 'StepTolerance', 1e-12);

%% Single shooting: unknowns are the two initial costates [lam_r0; lam_v0]
lam_guess = [0; 0];   % away from the analytic (3/2, 3/2)
[lam_sol, res, ef] = fsolve(@(L) sled_residual(L, tf, rf, vf, opts_ode), ...
                            lam_guess, opts_fs);

%% Report
lam_r0 = lam_sol(1); lam_v0 = lam_sol(2);
fprintf('=== Rocket sled shooting validation ===\n');
fprintf('  exitflag         : %d\n', ef);
fprintf('  recovered lam_r0 : %.10f   (analytic 1.5)\n', lam_r0);
fprintf('  recovered lam_v0 : %.10f   (analytic 1.5)\n', lam_v0);
fprintf('  terminal residual: ||[r(tf)-rf; v(tf)-vf]|| = %.3e\n', norm(res));
fprintf('  costate error    : ||lam - [1.5;1.5]||       = %.3e\n', norm(lam_sol - 1.5));

%% Cross-check the recovered control against the analytic u*(t) = (3/4)(1 - t)
[tt, S] = ode45(@(t,s) sled_ode(t,s), [0 tf], [0;0;lam_r0;lam_v0], opts_ode);
u_num = S(:,4)/2;                 % u = lam_v/2
u_ana = 0.75*(1 - tt);
fprintf('  max |u_num - u*| : %.3e\n', max(abs(u_num - u_ana)));

if ef > 0 && norm(lam_sol - 1.5) < 1e-6
    fprintf('\nPASS: shooting recovers the closed-form optimum.\n');
else
    fprintf('\nFAIL: solver did not recover the analytic optimum.\n');
end

%% ===================== LOCAL FUNCTIONS =====================
function ds = sled_ode(~, s)
% Sled RHS, u = lam_v/2.
%   INPUT
%     s - state [r; v; lam_r; lam_v]
%   OUTPUT
%     ds - derivative (4x1)
    ds = [ s(2);       % rdot = v
           s(4)/2;     % vdot = u = lam_v/2
           0;          % lam_r_dot = 0
          -s(3) ];     % lam_v_dot = -lam_r
end

function res = sled_residual(L, tf, rf, vf, opts_ode)
% Shooting residual at tf, integrating from rest.
%   INPUT
%     L        - initial costates [lam_r0; lam_v0]
%     tf       - final time
%     rf, vf   - targets at tf
%     opts_ode - ode45 options
%   OUTPUT
%     res - [r(tf)-rf; v(tf)-vf]
    [~, S] = ode45(@(t,s) sled_ode(t,s), [0 tf], [0; 0; L(1); L(2)], opts_ode);
    res = [S(end,1) - rf; S(end,2) - vf];
end
