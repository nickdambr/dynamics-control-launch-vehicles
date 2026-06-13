function dx = ode_lpv_flex(t, x, M)
%ODE_LPV_FLEX  LTV flexible pitch-plane RHS (bending + TVC + varying notch).
%
%   dx = ODE_LPV_FLEX(t, x, M) is the right-hand side of the time-varying
%   FLEXIBLE launch-vehicle model used by the T008 showcase: the HM3 Task-2
%   architecture (6-state plant with the first bending mode, INS coupling,
%   TVC actuator + transport delay, bending notch) made LPV over the ascent.
%   The 13-state vector stacks plant, notch and actuator states:
%
%       x = [ z zdot theta thetadot eta etadot | xn1 xn2 | x_tvc(1:nt) ]'
%
%   Plant (BUILD_PLANT_FULL form, all coefficients time-varying):
%       zddot     = a1*zdot + (a1*V+a4)*theta + a3*delta - a1*V*alpha_w
%       thetaddot = (A6/V)*zdot + A6*theta    + K1*delta - A6 *alpha_w
%       etaddot   = -omega^2*eta - 2*zBM*omega*etadot + aqk*delta
%   INS measurements (bending leaks in through sigma_ins(t), phi_ins(t)):
%       theta_m = theta + sigma*eta ,  z_m = z - phi*eta   (and the rates)
%   Controller (frozen max-q gains, or the schedule): u_pd from the meas{}.
%   Varying notch (controllable-canonical realisation, centre wn(t)):
%       xn1' = xn2 ;  xn2' = -wn^2*xn1 - 2*zD*wn*xn2 + u_pd
%       v    = 2*(zN - zD)*wn*xn2 + u_pd
%   TVC (LTI servo x Pade delay, state space M.tvc): delta = Ct*x_tvc + Dt*v.
%
%   M packs the INIT_SIMULINK_LPV interpolants (fa1..fK1, fomega, faqk, fsig,
%   fphi, windfun, fKp/fKd), the notch centre handle M.fwn (fixed = omega(72),
%   varying = omega(t)), the deep-notch constants (zN, zD, zBM), the TVC state
%   space (At,Bt,Ct,Dt), the frozen gains and the logical M.sched.
%
%   No input validation by design (ode45 inner loop). See ODE_LPV_ASCENT.
%
%   See also INIT_SIMULINK_LPV, MAIN_FLEX, BUILD_PLANT_FULL, BUILD_NOTCH_FILTER.

xp = x(1:6);  xn = x(7:8);  xt = x(9:end);
z = xp(1); zdot = xp(2); theta = xp(3); thetadot = xp(4); eta = xp(5); etadot = xp(6);

% --- time-varying coefficients + wind ---
a1 = M.fa1(t); a3 = M.fa3(t); a4 = M.fa4(t); A6 = M.fA6(t); K1 = M.fK1(t); V = M.fV(t);
w  = M.fomega(t); aqk = M.faqk(t); sig = M.fsig(t); phi = M.fphi(t);
aw = M.windfun(t);

% --- INS measurements (bending-contaminated) and PD command ---
theta_m = theta + sig*eta;  thetadot_m = thetadot + sig*etadot;
z_m     = z - phi*eta;      zdot_m     = zdot - phi*etadot;
if M.sched, Kp = M.fKp(t); Kd = M.fKd(t); else, Kp = M.Kp_th0; Kd = M.Kd_th0; end
u_pd = -(Kp*theta_m + Kd*thetadot_m + M.Kp_z*z_m + M.Kd_z*zdot_m);

% --- varying notch  ->  TVC  ->  physical deflection delta ---
wn = M.fwn(t);
v  = 2*(M.zN - M.zD)*wn*xn(2) + u_pd;
xn1d = xn(2);
xn2d = -wn^2*xn(1) - 2*M.zD*wn*xn(2) + u_pd;
delta = M.Ct*xt + M.Dt*v;
xtd   = M.At*xt + M.Bt*v;

% --- plant derivatives ---
zdd   = a1*zdot + (a1*V + a4)*theta + a3*delta - a1*V*aw;
thdd  = (A6/V)*zdot + A6*theta + K1*delta - A6*aw;
etadd = -w^2*eta - 2*M.zBM*w*etadot + aqk*delta;

dx = [zdot; zdd; thetadot; thdd; etadot; etadd; xn1d; xn2d; xtd];
end
