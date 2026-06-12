function dz = ode_burn(t, z, p)
% ODE for powered flight with linear tangent steering law
%   z = [x; y; vx; vy; m; lam_m]
%   p = struct with fields: T, Q, c, lam_vx0, lam_vy0, lam_y
%
% Costate structure (linear tangent law):
%   lam_x  = 0        (x free at tf)
%   lam_y  = const    = p.lam_y
%   lam_vx = const    = p.lam_vx0
%   lam_vy = linear   = p.lam_vy0 - p.lam_y * t
%
% Optimal thrust angle: phi = atan2(lam_vy, lam_vx)
%
% No arguments validation by design: this RHS is called ~1e6+ times by
% ode45 inside the fsolve shooting loops; validate at the call site.

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
