function dz = ode_burn(t, z, p)
% Powered-flight RHS, linear-tangent steering.
%   INPUT
%     t - time (scalar)
%     z - state [x; y; vx; vy; m; lam_m]
%     p - struct: T, Q, c, lam_vx0, lam_vy0, lam_y
%   OUTPUT
%     dz - state derivative (6x1)
%
% Costates: lam_x=0, lam_y=p.lam_y, lam_vx=p.lam_vx0,
% lam_vy=p.lam_vy0 - p.lam_y*t. Thrust angle phi=atan2(lam_vy,lam_vx).
%
% No arguments block by design: called ~1e6+ times inside ode45/fsolve;
% validate at the call site.

vx = z(3); vy = z(4); m = z(5);

% Costates
lam_vx = p.lam_vx0;
lam_vy = p.lam_vy0 - p.lam_y * t;
lam_v_norm = sqrt(lam_vx^2 + lam_vy^2);

% Optimal thrust angle
phi = atan2(lam_vy, lam_vx);

% State derivatives
dz = zeros(6,1);
dz(1) = vx;                          % dx/dt
dz(2) = vy;                          % dy/dt
dz(3) = (p.T / m) * cos(phi);        % dvx/dt
dz(4) = (p.T / m) * sin(phi) - 1;    % dvy/dt  (g = 1 nondim)
dz(5) = -p.Q;                         % dm/dt
dz(6) = (p.T / m^2) * lam_v_norm;    % dlam_m/dt

end
