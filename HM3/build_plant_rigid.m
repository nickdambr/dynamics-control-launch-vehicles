function G = build_plant_rigid(p)
%BUILD_PLANT_RIGID  Rigid-body pitch-plane LV plant at max-qbar (Task 1).
%
%   G = BUILD_PLANT_RIGID(p) returns a state-space model of the rigid LV
%   obtained from Eq. (1) of the assignment by dropping the bending mode
%   (eta, etadot). The 4 states are
%
%       x = [z, zdot, theta, thetadot]'
%
%   with inputs u = [delta, alpha_w]' (TVC deflection and wind angle of
%   attack) and outputs
%
%       y = [theta_m, thetadot_m, z_m, zdot_m, theta, z, zdot]'
%
%   The first four are the INS measurements fed back to the controller
%   (identical to the true states in the rigid case, since there is no
%   bending contamination) and the last three are convenience signals for
%   plotting. Parameters come from LOAD_HW3_PARAMS.
%
%   See also BUILD_PLANT_FULL, LOAD_HW3_PARAMS.

A = [0   1        0          0;
     0   p.a1     p.a1*p.V+p.a4   0;
     0   0        0          1;
     0   p.A6/p.V p.A6       0];

Bd = [0; p.a3; 0; p.K1];          % delta column
Bw = [0; -p.a1*p.V; 0; -p.A6];    % alpha_w column

% Outputs: [theta_m, thetadot_m, z_m, zdot_m, theta, z, zdot]
%          state index:   z=1  zdot=2  theta=3  thetadot=4
C = [0 0 1 0;    % theta_m   = theta
     0 0 0 1;    % thetadot_m= thetadot
     1 0 0 0;    % z_m       = z
     0 1 0 0;    % zdot_m    = zdot
     0 0 1 0;    % theta
     1 0 0 0;    % z
     0 1 0 0];   % zdot
D = zeros(7, 2);

G = ss(A, [Bd Bw], C, D);
G.StateName  = {'z','zdot','theta','thetadot'};
G.InputName  = {'delta','alpha_w'};
G.OutputName = {'theta_m','thetadot_m','z_m','zdot_m','theta','z','zdot'};
G.Name = 'Plant_Rigid';
end
