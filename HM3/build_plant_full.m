function G = build_plant_full(p, meas)
% Full pitch-plane LV plant with bending mode (Task 2/3, Eq. 1).
%   INPUT
%     p    - param struct (load_hw3_params)
%     meas - 'ins' (Eq. 2, bending leaks into measurements) | 'true'
%            (uncontaminated feedback). Default 'ins'.
%   OUTPUT
%     G - ss, 6 states [z zdot theta thetadot eta etadot], in [delta alpha_w],
%         out [theta_m thetadot_m z_m zdot_m theta z zdot]
%   INS bending contamination (sigma_ins, phi_ins) is what destabilises the
%   loop at wBM and motivates the Task-2 notch.

arguments
    p (1,1) struct
    meas {mustBeTextScalar} = 'ins'   % 'ins' | 'true', case-insensitive
end

w = p.wBM;  z = p.zBM;

A = [0   1            0              0   0      0;
     0   p.a1         p.a1*p.V+p.a4  0   0      0;
     0   0            0              1   0      0;
     0   p.A6/p.V     p.A6           0   0      0;
     0   0            0              0   0      1;
     0   0            0              0  -w^2   -2*z*w];

Bd = [0; p.a3; 0; p.K1; 0; -p.phi_tvc*p.Tc];   % delta column
Bw = [0; -p.a1*p.V; 0; -p.A6; 0; 0];           % alpha_w column

% state index:  z=1 zdot=2 theta=3 thetadot=4 eta=5 etadot=6
switch lower(meas)
    case 'ins'   % Eq. (2): bending leaks into the measurements
        Cm = [0 0 1 0  p.sigma_ins 0;          % theta_m   = theta + sigma*eta
              0 0 0 1  0 p.sigma_ins;          % thetadot_m= thetadot + sigma*etadot
              1 0 0 0 -p.phi_ins 0;            % z_m       = z - phi*eta
              0 1 0 0  0 -p.phi_ins];          % zdot_m    = zdot - phi*etadot
    case 'true'  % uncontaminated feedback
        Cm = [0 0 1 0 0 0;
              0 0 0 1 0 0;
              1 0 0 0 0 0;
              0 1 0 0 0 0];
    otherwise
        error('build_plant_full:meas', 'meas must be ''ins'' or ''true''.');
end
Cplot = [0 0 1 0 0 0;     % theta (true)
         1 0 0 0 0 0;     % z (true)
         0 1 0 0 0 0];    % zdot (true)
C = [Cm; Cplot];
D = zeros(7, 2);

G = ss(A, [Bd Bw], C, D);
G.StateName  = {'z','zdot','theta','thetadot','eta','etadot'};
G.InputName  = {'delta','alpha_w'};
G.OutputName = {'theta_m','thetadot_m','z_m','zdot_m','theta','z','zdot'};
G.Name = 'Plant_Full';
end
