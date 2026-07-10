function r = simulate_gust_response(T, w)
% Closed-loop time response to a wind gust (theta_ref = 0).
%   INPUT
%     T - closed loop (assemble_loop), in {alpha_w theta_ref}, out {theta z zdot delta}
%     w - wind struct (load_wind_profile)
%   OUTPUT
%     r - struct: t, alphaw, theta, z, zdot, delta, alpha = theta+zdot/V+alphaw
%         (aero-load driver at max-qbar), plus peak_* metrics

arguments
    T {mustBeA(T, 'lti')}
    w (1,1) struct
end

t = w.t(:);
u = [w.alphaw(:), zeros(numel(t),1)];   % [alpha_w, theta_ref]
y = lsim(T, u, t);

r.t      = t;
r.alphaw = w.alphaw(:);
r.theta  = y(:,1);
r.z      = y(:,2);
r.zdot   = y(:,3);
r.delta  = y(:,4);
r.alpha  = r.theta + r.zdot/w.V + r.alphaw;   % total angle of attack

r.peak_theta = max(abs(r.theta));
r.peak_z     = max(abs(r.z));
r.peak_delta = max(abs(r.delta));
r.peak_alpha = max(abs(r.alpha));
end
