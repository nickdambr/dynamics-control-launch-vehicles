function [L, T, info] = assemble_loop(G, K, Wact)
%ASSEMBLE_LOOP  Close the PD attitude/drift loop around the LV plant.
%
%   [L, T] = ASSEMBLE_LOOP(G, K) builds the negative-feedback control loop
%
%       u_pd  = Kp_th*(theta_ref - theta_m) - Kd_th*thetadot_m
%               - Kp_z*z_m - Kd_z*zdot_m
%       delta = Wact(s) * u_pd
%
%   around the plant G (from BUILD_PLANT_RIGID / BUILD_PLANT_FULL). The
%   gains live in struct K with fields Kp_th, Kd_th, Kp_z, Kd_z (the lateral
%   gains Kp_z, Kd_z are expected to be small and negative, per the
%   assignment guidelines).
%
%   Outputs:
%     L  - SISO open-loop transfer broken at the TVC deflection 'delta',
%          in the 1+L convention (use MARGIN(L) / NICHOLS(L)).
%     T  - closed-loop model from inputs {alpha_w, theta_ref} to outputs
%          {theta, z, zdot, delta}, for time-domain simulation.
%
%   [L, T] = ASSEMBLE_LOOP(G, K, Wact) inserts the LTI block Wact (e.g.
%   TVC + delay + notch, in series) between the controller and the plant.
%   Default Wact = 1 (rigid Task-1 case, ideal actuator).
%
%   See also BUILD_PLANT_RIGID, BUILD_PLANT_FULL, BUILD_TVC,
%   BUILD_NOTCH_FILTER, DESIGN_CONTROLLER.

arguments
    G {mustBeA(G, 'lti')}
    K (1,1) struct
    Wact = tf(1)
end

if isempty(Wact), Wact = tf(1); end    % [] is the documented ideal-actuator alias

% --- controller as a static gain block (named IO) ---
% u_pd = [Kp_th, -Kp_th, -Kd_th, -Kp_z, -Kd_z] * [theta_ref;theta_m;thetadot_m;z_m;zdot_m]
Kc = ss([K.Kp_th, -K.Kp_th, -K.Kd_th, -K.Kp_z, -K.Kd_z]);
Kc.InputName  = {'theta_ref','theta_m','thetadot_m','z_m','zdot_m'};
Kc.OutputName = {'u_pd'};

% --- actuator / filter chain u_pd -> delta ---
Wa = ss(Wact);
Wa.InputName  = {'u_pd'};
Wa.OutputName = {'delta'};

% --- closed loop, retaining 'delta' as an analysis point ---
T = connect(G, Kc, Wa, {'alpha_w','theta_ref'}, {'theta','z','zdot','delta'}, {'delta'});

% --- open-loop transfer at the break point 'delta' (1+L convention) ---
L = getLoopTransfer(T, 'delta', -1);
L = minreal(tf(L), 1e-6);

info = struct('Kc', Kc, 'Wact', Wa);
end
