function dx = ode_descent(x, u, Vc)
% Non-dim powered-descent point-mass RHS. Gravity = -1.
%   INPUT
%     x  - state [x; y; vx; vy; m]
%     u  - thrust [Tx; Ty]
%     Vc - V_ref/c (Tsiolkovsky number)
%   OUTPUT
%     dx - d/dt [x; y; vx; vy; m] (5x1)
%
% No arguments validation by design: hot-loop RHS inside ode45/fmincon;
% validate at the call site.

Tmag = sqrt(u(1)^2 + u(2)^2);
dx = [ x(3); x(4); u(1)/x(5); u(2)/x(5) - 1; -Vc * Tmag ];

end
