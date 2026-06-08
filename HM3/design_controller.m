function [K, m] = design_controller(G, Wact, varargin)
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

ip = inputParser;
ip.addParameter('Kp_z',-1e-3);
ip.addParameter('Kd_z',-1e-3);
ip.addParameter('GM',6);
ip.addParameter('PM',30);
ip.addParameter('K0',[2.0 1.4]);
ip.addParameter('verbose',true);
ip.parse(varargin{:});
o = ip.Results;

if nargin < 2, Wact = tf(1); end
if isempty(Wact), Wact = tf(1); end

    function c = cost(x)
        Kt.Kp_th = exp(x(1));  Kt.Kd_th = exp(x(2));
        Kt.Kp_z  = o.Kp_z;     Kt.Kd_z  = o.Kd_z;
        warning('off','Control:analysis:MarginUnstable');
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

warning('off','Control:analysis:MarginUnstable');
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
