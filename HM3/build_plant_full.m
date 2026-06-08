function G = build_plant_full(p, meas)
%BUILD_PLANT_FULL  Full pitch-plane LV plant with bending mode (Task 2/3).
%
%   G = BUILD_PLANT_FULL(p) returns the 6-state model of Eq. (1):
%
%       x = [z, zdot, theta, thetadot, eta, etadot]'
%
%   with inputs u = [delta, alpha_w]' and outputs
%
%       y = [theta_m, thetadot_m, z_m, zdot_m, theta, z, zdot]'
%
%   The measurement block (first four outputs) follows the INS model of
%   Eq. (2): the bending generalised coordinate eta leaks into the gyro and
%   accelerometer channels through sigma_ins and phi_ins. It is precisely
%   this contamination that destabilises the loop at the bending frequency
%   and motivates the notch / lead-lag filter of Task 2.
%
%   G = BUILD_PLANT_FULL(p, 'true') bypasses the INS model and feeds back
%   the true (uncontaminated) states instead, useful as a debugging
%   baseline. Default is 'ins'.
%
%   See also BUILD_PLANT_RIGID, BUILD_INS_MODEL, BUILD_NOTCH_FILTER.

if nargin < 2 || isempty(meas), meas = 'ins'; end

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
