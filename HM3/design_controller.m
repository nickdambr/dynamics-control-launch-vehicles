function [K, m] = design_controller(G, Wact, o)
% Design the pitch PD gains on the FULL launch-vehicle loop (D'Antuono method).
%
%   Starting from the canonical closed-form gains on the decoupled rotational
%   dynamics (D'Antuono Eq. 3.6-3.7) Kp = 2*A6/K1, Kd = sqrt(A6)/K1, the gains
%   are RE-TUNED on the FULL loop (attitude + lateral drift + actuator) so the
%   CLASSIFIED Aero GM and Rigid PM hit the assignment targets. This retune is
%   what actually meets the requirement: the lateral-drift feedback erodes the
%   aerodynamic gain margin (the canonical decoupled 6 dB drops to ~4 dB on the
%   full loop). Margins are classified by frequency band (classify_margins);
%   the closed-loop stability verdict is isstable() of the full loop.
%
%   INPUT
%     G    - rigid plant (build_plant_rigid); A6,K1 read by state name
%     Wact - actuator chain; [] or tf(1) = ideal actuator (Task 1)
%     o    - name-value:
%              Kp_z, Kd_z  fixed lateral-drift gains (default -1e-3, -1e-3)
%              GM, PM      margin targets            (default 6 dB, 30 deg)
%              K0          ignored (kept for call compatibility)
%              verbose     print result              (default true)
%   OUTPUT
%     K - tuned gain struct (Kp_th, Kd_th, Kp_z, Kd_z)
%     m - classified margins (see classify_margins) + stable, L, T

arguments
    G {mustBeA(G, 'lti')}
    Wact = tf(1)
    o.Kp_z    (1,1) {mustBeNumeric, mustBeReal} = -1e-3
    o.Kd_z    (1,1) {mustBeNumeric, mustBeReal} = -1e-3
    o.GM      (1,1) {mustBeNumeric, mustBeReal, mustBePositive} = 6
    o.PM      (1,1) {mustBeNumeric, mustBeReal, mustBePositive} = 30
    o.K0      (1,2) {mustBeNumeric, mustBeReal} = [0 0]   % accepted, unused
    o.w_flex    (1,1) {mustBeNumeric, mustBeReal, mustBePositive} = Inf  % rigid/flex bound
    o.w_flex_hi (1,1) {mustBeNumeric, mustBeReal, mustBePositive} = Inf  % upper flex bound
    o.w_bending (1,1) {mustBeNumeric, mustBeReal} = NaN                  % bending freq
    o.verbose (1,1) logical = true
end
if isempty(Wact), Wact = tf(1); end

% decoupled rotational coefficients (by state/input name)
iTh = strcmp(G.StateName, 'theta');
iTd = strcmp(G.StateName, 'thetadot');
iDe = strcmp(G.InputName, 'delta');
A6  = G.A(iTd, iTh);
K1  = G.B(iTd, iDe);
w_drift = 0.3*sqrt(A6);          % drift/rigid boundary for the classifier
bands = {'w_drift', w_drift, 'w_flex', o.w_flex, 'w_flex_hi', o.w_flex_hi, ...
         'w_bending', o.w_bending};   % Task 1: defaults (no bending); Task 2: full-loop bands

% MARGIN warns on every conditionally-stable evaluation; mute for the search.
warnState   = warning('off', 'Control:analysis:MarginUnstable');
restoreWarn = onCleanup(@() warning(warnState));

% Re-tune (Kp,Kd) on the full loop from the canonical closed-form start.
x0 = log([2*A6/K1, sqrt(A6)/K1]);         % D'Antuono Eq. 3.6-3.7
xo = fminsearch(@cost, x0, ...
                optimset('Display','off','TolX',1e-4,'TolFun',1e-3,'MaxFunEvals',400));

K.Kp_th = exp(xo(1));
K.Kd_th = exp(xo(2));
K.Kp_z  = o.Kp_z;
K.Kd_z  = o.Kd_z;

[L, T]   = assemble_loop(G, K, Wact);
L        = minreal(L, 1e-6);
m        = classify_margins(L, bands{:});
m.stable = isstable(T);
m.L      = L;
m.T      = T;

if o.verbose
    fprintf(['  PD design (full loop): Kp_th=%.4f Kd_th=%.4f | Kp_z=%.1e Kd_z=%.1e\n' ...
             '    Aero |GM|=%.2f dB @%.2f rad/s  Rigid PM=%.1f deg @%.2f rad/s  DM=%.0f ms | CL stable: %d\n'], ...
            K.Kp_th, K.Kd_th, K.Kp_z, K.Kd_z, ...
            abs(m.aeroGM_dB), m.aeroGM_w, m.rigidPM_deg, m.rigidPM_w, 1e3*m.DM_s, m.stable);
end

    function c = cost(x)
        % Full-loop margin-matching cost. x = log([Kp_th Kd_th]).
        Kt.Kp_th = exp(x(1));  Kt.Kd_th = exp(x(2));
        Kt.Kp_z  = o.Kp_z;     Kt.Kd_z  = o.Kd_z;
        [Lt, Tt] = assemble_loop(G, Kt, Wact);
        Lt = minreal(Lt, 1e-6);
        mt = classify_margins(Lt, bands{:});
        if isnan(mt.aeroGM_dB) || isnan(mt.rigidPM_deg)
            c = 1e6;  return;                 % lost a required crossover
        end
        c = (abs(mt.aeroGM_dB) - o.GM)^2 + (mt.rigidPM_deg - o.PM)^2;
        if ~isstable(Tt), c = c + 1e4; end
    end
end
