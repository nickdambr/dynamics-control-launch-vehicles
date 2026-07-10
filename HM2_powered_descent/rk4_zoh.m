function x_next = rk4_zoh(x, u, dt, Vc, n_sub)
% RK4 propagation of ode_descent over one ZOH interval, u held constant.
%   INPUT
%     x     - state [x; y; vx; vy; m]
%     u     - thrust [Tx; Ty] (held over dt)
%     dt    - interval length
%     Vc    - V_ref/c
%     n_sub - RK4 substeps
%   OUTPUT
%     x_next - state at end of interval (5x1)
%
% No arguments validation by design: hot-loop propagator inside fmincon;
% validate at the call site.

h = dt / n_sub;
for ii = 1:n_sub
    k1 = ode_descent(x,            u, Vc);
    k2 = ode_descent(x + 0.5*h*k1, u, Vc);
    k3 = ode_descent(x + 0.5*h*k2, u, Vc);
    k4 = ode_descent(x +     h*k3, u, Vc);
    x  = x + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
end
x_next = x;

end
