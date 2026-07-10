function dx = ode_descent_uacc(x, uacc, Vc)
% Non-dim powered-descent RHS with the ACCELERATION u = T/m held constant
% (the ZOH convention native to the GFOLD log-mass transcription). Gravity = -1.
%   INPUT
%     x    - state [x; y; vx; vy; m]
%     uacc - acceleration [ux; uy] = T/m (held constant over the interval)
%     Vc   - V_ref/c (Tsiolkovsky number)
%   OUTPUT
%     dx   - d/dt [x; y; vx; vy; m] (5x1)
%
% With u held constant the thrust T = m(t)*u floats with the (depleting)
% mass, so vx_dot = ux and vy_dot = uy - 1 are exact, while the mass row
% reads m_dot = -Vc*||T|| = -Vc*m*||u||. Contrast ode_descent.m, which holds
% the thrust vector T (not the acceleration) piecewise constant.
%
% No arguments validation by design: hot-loop RHS inside ode45.

umag = sqrt(uacc(1)^2 + uacc(2)^2);
dx = [ x(3); x(4); uacc(1); uacc(2) - 1; -Vc * x(5) * umag ];

end
