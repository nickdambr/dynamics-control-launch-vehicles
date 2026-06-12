function x_next = rk4_zoh(x, u, dt, Vc, n_sub)
% Propagate the non-dim descent dynamics (ode_descent) over one ZOH
% interval of length dt with a fixed-step RK4 scheme, holding the control
% u constant, using n_sub substeps.
% Shared by main_task2.m (nonlinear ZOH transcription) and the test suite.
%
% No arguments validation by design: hot-loop propagator called by the
% fmincon ZOH transcription; validate at the call site.

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
