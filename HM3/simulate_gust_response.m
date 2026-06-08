function r = simulate_gust_response(T, w)
%SIMULATE_GUST_RESPONSE  Closed-loop time response to a wind gust.
%
%   r = SIMULATE_GUST_RESPONSE(T, w) simulates the closed-loop model T
%   (from ASSEMBLE_LOOP, inputs {alpha_w, theta_ref}, outputs
%   {theta, z, zdot, delta}) driven by the wind disturbance w (from
%   LOAD_WIND_PROFILE), with theta_ref = 0.
%
%   The output struct r contains the time vector and the key time
%   histories requested by the assignment: theta, z, zdot, delta, plus the
%   driving alpha_w and scalar peak metrics.
%
%   See also LOAD_WIND_PROFILE, ASSEMBLE_LOOP.

t = w.t(:);
u = [w.alphaw(:), zeros(numel(t),1)];   % [alpha_w, theta_ref]
y = lsim(T, u, t);

r.t      = t;
r.alphaw = w.alphaw(:);
r.theta  = y(:,1);
r.z      = y(:,2);
r.zdot   = y(:,3);
r.delta  = y(:,4);

r.peak_theta = max(abs(r.theta));
r.peak_z     = max(abs(r.z));
r.peak_delta = max(abs(r.delta));
end
