function G = build_plant_rigid(p)
% Rigid pitch-plane LV plant at max-qbar (Task 1, Eq. 1 minus bending).
%   INPUT
%     p - param struct (load_hw3_params)
%   OUTPUT
%     G - ss, 4 states [z zdot theta thetadot], in [delta alpha_w],
%         out [theta_m thetadot_m z_m zdot_m theta z zdot]
%   First 4 outputs feed the controller (= true states here, no bending);
%   last 3 are plotting signals.

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
