function [L, T, info] = assemble_loop(G, K, Wact)
% Close the PD attitude/drift loop around the LV plant.
%   u_pd  = Kp_th*(theta_ref-theta_m) - Kd_th*thetadot_m - Kp_z*z_m - Kd_z*zdot_m
%   delta = Wact(s) * u_pd
%   INPUT
%     G    - plant (build_plant_rigid / build_plant_full)
%     K    - gain struct: Kp_th, Kd_th, Kp_z, Kd_z (lateral gains small/negative)
%     Wact - series block controller->plant (TVC+delay+notch). Default 1
%            (ideal actuator); [] is the same alias.
%   OUTPUT
%     L    - SISO open loop broken at 'delta', 1+L convention (margin/nichols)
%     T    - closed loop in {alpha_w, theta_ref} -> out {theta z zdot delta}
%     info - struct with the controller/actuator IO blocks

arguments
    G {mustBeA(G, 'lti')}
    K (1,1) struct
    Wact = tf(1)
end

if isempty(Wact), Wact = tf(1); end    % [] => ideal actuator

% --- controller as a static gain block (named IO) ---
% u_pd = [Kp_th -Kp_th -Kd_th -Kp_z -Kd_z] * [theta_ref;theta_m;thetadot_m;z_m;zdot_m]
Kc = ss([K.Kp_th, -K.Kp_th, -K.Kd_th, -K.Kp_z, -K.Kd_z]);
Kc.InputName  = {'theta_ref','theta_m','thetadot_m','z_m','zdot_m'};
Kc.OutputName = {'u_pd'};

% --- actuator / filter chain u_pd -> delta ---
Wa = ss(Wact);
Wa.InputName  = {'u_pd'};
Wa.OutputName = {'delta'};

% --- closed loop, keeping 'delta' as an analysis point ---
T = connect(G, Kc, Wa, {'alpha_w','theta_ref'}, {'theta','z','zdot','delta'}, {'delta'});

% --- open loop at the break point 'delta' (1+L convention) ---
L = getLoopTransfer(T, 'delta', -1);
L = minreal(tf(L), 1e-6);

info = struct('Kc', Kc, 'Wact', Wa);
end
