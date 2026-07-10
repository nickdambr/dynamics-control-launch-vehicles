function dx = ode_lpv_ascent(t, x, M)
% LTV rigid pitch-plane RHS, full-ascent LPV baseline. Ideal-actuator PD
% closed inside with theta_ref = 0. Coefficients c1..c7 are time interpolants.
%   INPUT
%     t - flight time [s]
%     x - state [z; zdot; theta; thetadot]
%     M - struct: fc1..fc7, windfun, fKp/fKd, frozen gains, sched (0/1)
%   OUTPUT
%     dx - state derivative (4x1)
%
% No arguments validation by design: ode45 inner loop (~1e5+ calls).
% Validation lives in INIT_SIMULINK_LPV / MAIN_FULL_ASCENT.

aw = M.windfun(t);                       % wind angle of attack alpha_w(t)

% controller: frozen max-qbar gains, or the schedule
if M.sched
    Kp = M.fKp(t);  Kd = M.fKd(t);
else
    Kp = M.Kp_th0;  Kd = M.Kd_th0;
end
delta = -(Kp*x(3) + Kd*x(4) + M.Kp_z*x(1) + M.Kd_z*x(2));

% time-varying rigid plant
dx = [ x(2);
       M.fc1(t)*x(2) + M.fc2(t)*x(3) + M.fc3(t)*delta - M.fc4(t)*aw;
       x(4);
       M.fc5(t)*x(2) + M.fc6(t)*x(3) + M.fc7(t)*delta - M.fc6(t)*aw ];
end
