function dx = ode_lpv_flex(t, x, M)
% LTV flexible pitch-plane RHS: bending + INS coupling + TVC/delay + varying
% notch (notch centre M.fwn: omega(72) fixed, omega(t) varying). All plant
% coefficients are time interpolants. PD frozen at max-q (M.sched typ. false).
%   INPUT
%     t - flight time [s]
%     x - state [z zdot theta thetadot eta etadot | xn1 xn2 | x_tvc(1:nt)]
%     M - struct: fa1..fK1, fV, fomega, faqk, fsig, fphi, windfun, fwn,
%                 zN/zD/zBM, At/Bt/Ct/Dt (TVC ss), fKp/fKd, frozen gains, sched
%   OUTPUT
%     dx - state derivative (13x1, nt = TVC order)
%
% No arguments validation by design: ode45 inner loop. See ODE_LPV_ASCENT.

xp = x(1:6);  xn = x(7:8);  xt = x(9:end);
z = xp(1); zdot = xp(2); theta = xp(3); thetadot = xp(4); eta = xp(5); etadot = xp(6);

% time-varying coefficients + wind
a1 = M.fa1(t); a3 = M.fa3(t); a4 = M.fa4(t); A6 = M.fA6(t); K1 = M.fK1(t); V = M.fV(t);
w  = M.fomega(t); aqk = M.faqk(t); sig = M.fsig(t); phi = M.fphi(t);
aw = M.windfun(t);

% INS measurements (bending-contaminated) and PD command
theta_m = theta + sig*eta;  thetadot_m = thetadot + sig*etadot;
z_m     = z - phi*eta;      zdot_m     = zdot - phi*etadot;
if M.sched, Kp = M.fKp(t); Kd = M.fKd(t); else, Kp = M.Kp_th0; Kd = M.Kd_th0; end
u_pd = -(Kp*theta_m + Kd*thetadot_m + M.Kp_z*z_m + M.Kd_z*zdot_m);

% varying notch -> TVC -> physical deflection delta
wn = M.fwn(t);
v  = 2*(M.zN - M.zD)*wn*xn(2) + u_pd;
xn1d = xn(2);
xn2d = -wn^2*xn(1) - 2*M.zD*wn*xn(2) + u_pd;
delta = M.Ct*xt + M.Dt*v;
xtd   = M.At*xt + M.Bt*v;

% plant derivatives
zdd   = a1*zdot + (a1*V + a4)*theta + a3*delta - a1*V*aw;
thdd  = (A6/V)*zdot + A6*theta + K1*delta - A6*aw;
etadd = -w^2*eta - 2*M.zBM*w*etadot + aqk*delta;

dx = [zdot; zdd; thetadot; thdd; etadot; etadd; xn1d; xn2d; xtd];
end
