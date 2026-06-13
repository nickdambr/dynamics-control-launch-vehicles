function [K, m] = design_controller(G, Wact, o)
%DESIGN_CONTROLLER  Tune the pitch PD gains to target gain/phase margins.
%
%   [K, m] = DESIGN_CONTROLLER(G, Wact) tunes the pitch proportional and
%   derivative gains (Kp_th, Kd_th) so that the open-loop transfer of
%   ASSEMBLE_LOOP attains the assignment targets
%
%       |GM| ~ 6 dB ,   |PM| ~ 30 deg
%
%   while keeping the closed loop stable. Because the rigid airframe is
%   open-loop unstable (pole at +sqrt(A6)), the loop is conditionally
%   stable and MARGIN reports the low-frequency aerodynamic gain margin;
%   the targets are therefore matched in magnitude (as read off the Nichols
%   chart). The lateral-drift gains Kp_z, Kd_z are small and negative per
%   the assignment guidelines and are held fixed during the search.
%
%   Pass Wact = [] (or tf(1)) for the ideal-actuator rigid case (Task 1),
%   or the TVC+delay(+notch) chain for the full model (Task 2).
%
%   Name/value options:
%     'Kp_z','Kd_z'   fixed lateral gains      (default -1e-3, -1e-3)
%     'GM','PM'       margin targets           (default 6 dB, 30 deg)
%     'K0'            [Kp_th0 Kd_th0] guess    (default [2.0 1.4])
%     'verbose'       print result             (default true)
%
%   See also ASSEMBLE_LOOP, MARGIN.

arguments
    G {mustBeA(G, 'lti')}
    Wact = tf(1)
    o.Kp_z    (1,1) {mustBeNumeric, mustBeReal} = -1e-3
    o.Kd_z    (1,1) {mustBeNumeric, mustBeReal} = -1e-3
    o.GM      (1,1) {mustBeNumeric, mustBeReal, mustBePositive} = 6
    o.PM      (1,1) {mustBeNumeric, mustBeReal, mustBePositive} = 30
    o.K0      (1,2) {mustBeNumeric, mustBeReal, mustBePositive} = [2.0 1.4]
    o.verbose (1,1) logical = true
end

if isempty(Wact), Wact = tf(1); end    % [] is the documented ideal-actuator alias

% MARGIN warns on every evaluation of the conditionally stable loop; mute it
% for the whole search and restore the caller's warning state on exit.
warnState = warning('off','Control:analysis:MarginUnstable');
restoreWarn = onCleanup(@() warning(warnState));

    function c = cost(x)
        Kt.Kp_th = exp(x(1));  Kt.Kd_th = exp(x(2));
        Kt.Kp_z  = o.Kp_z;     Kt.Kd_z  = o.Kd_z;
        [L,T] = assemble_loop(G, Kt, Wact);
        [Gm,Pm] = margin(L);
        gm_db = 20*log10(Gm);
        c = (abs(gm_db)-o.GM)^2 + (abs(Pm)-o.PM)^2;
        if ~isstable(T), c = c + 1e4; end           % keep CL stable
        if ~isfinite(c), c = 1e6; end
    end

x0  = log(o.K0);
opts = optimset('Display','off','TolX',1e-4,'TolFun',1e-4,'MaxFunEvals',400);
xopt = fminsearch(@cost, x0, opts);

K.Kp_th = exp(xopt(1));
K.Kd_th = exp(xopt(2));
K.Kp_z  = o.Kp_z;
K.Kd_z  = o.Kd_z;

[L,T] = assemble_loop(G, K, Wact);
[Gm,Pm,Wcg,Wcp] = margin(L);
m = struct('GM_dB',20*log10(Gm),'PM_deg',Pm,'wc_gain',Wcg,'wc_phase',Wcp, ...
           'stable',isstable(T),'L',L,'T',T);

if o.verbose
    fprintf(['  PD design: Kp_th=%.4f Kd_th=%.4f | Kp_z=%.1e Kd_z=%.1e\n' ...
             '             GM=%.2f dB (w=%.3g) PM=%.1f deg (w=%.3g) | CL stable: %d\n'], ...
            K.Kp_th,K.Kd_th,K.Kp_z,K.Kd_z, m.GM_dB,m.wc_gain,m.PM_deg,m.wc_phase,m.stable);
end
end
