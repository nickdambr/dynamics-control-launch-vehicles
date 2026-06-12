function dx = ode_descent(x, u, Vc)
% Non-dimensional powered-descent point-mass dynamics (HM2).
%   x  = [x; y; vx; vy; m]   non-dim state
%   u  = [Tx; Ty]            non-dim thrust components
%   Vc = V_ref / c           effective Tsiolkovsky number
%
% Returns dx = d/dt [x; y; vx; vy; m] with non-dim gravity = -1.
% Shared by main_task1.m, main_task2.m and the test suite.
%
% No arguments validation by design: hot-loop RHS called by ode45 and the
% fmincon transcriptions; validate at the call site.

Tmag = sqrt(u(1)^2 + u(2)^2);
dx = [ x(3); x(4); u(1)/x(5); u(2)/x(5) - 1; -Vc * Tmag ];

end
