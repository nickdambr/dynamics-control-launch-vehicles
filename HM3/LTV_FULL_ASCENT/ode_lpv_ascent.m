function dx = ode_lpv_ascent(t, x, M)
%ODE_LPV_ASCENT  LTV rigid pitch-plane RHS for the full-ascent LPV baseline.
%
%   dx = ODE_LPV_ASCENT(t, x, M) is the right-hand side of the time-varying
%   rigid launch-vehicle model, the source-of-truth replay of the Simulink
%   model hm3_full_ascent.slx. The state is the same as BUILD_PLANT_RIGID,
%
%       x = [z, zdot, theta, thetadot]'
%
%   and the dynamics use the "effective coefficients" precomputed by
%   INIT_SIMULINK_LPV (each is one interpolant of flight time):
%
%       zddot     = c1*zdot + c2*theta + c3*delta - c4*alpha_w
%       thetaddot = c5*zdot + c6*theta + c7*delta - c6*alpha_w
%
%   The ideal-actuator PD law (delta = u_pd, Task-1 style) is closed inside
%   the RHS with theta_ref = 0:
%
%       delta = -(Kp_th*theta + Kd_th*thetadot + Kp_z*z + Kd_z*zdot)
%
%   M is the model struct from INIT_SIMULINK_LPV (griddedInterpolant handles
%   fc1..fc7, windfun, fKp/fKd, plus the scalar frozen gains and the logical
%   M.sched selecting frozen (0) vs gain-scheduled (1) pitch gains).
%
%   No input validation by design: this is an ode45 inner-loop function
%   (~1e5+ calls). Validation lives in INIT_SIMULINK_LPV / MAIN_FULL_ASCENT.
%
%   See also INIT_SIMULINK_LPV, MAIN_FULL_ASCENT, BUILD_PLANT_RIGID.

aw = M.windfun(t);                       % wind angle of attack alpha_w(t)

% --- controller (frozen max-qbar gains, or the time-varying schedule) ---
if M.sched
    Kp = M.fKp(t);  Kd = M.fKd(t);
else
    Kp = M.Kp_th0;  Kd = M.Kd_th0;
end
delta = -(Kp*x(3) + Kd*x(4) + M.Kp_z*x(1) + M.Kd_z*x(2));

% --- time-varying rigid plant ---
dx = [ x(2);
       M.fc1(t)*x(2) + M.fc2(t)*x(3) + M.fc3(t)*delta - M.fc4(t)*aw;
       x(4);
       M.fc5(t)*x(2) + M.fc6(t)*x(3) + M.fc7(t)*delta - M.fc6(t)*aw ];
end
