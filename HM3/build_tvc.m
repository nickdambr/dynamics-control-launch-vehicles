function Wtvc = build_tvc(p, padeOrder)
%BUILD_TVC  TVC actuator transfer function with pure delay (Eq. 3).
%
%   Wtvc = BUILD_TVC(p) returns the second-order thrust-vector-control
%   actuator model
%
%       W_TVC(s) = wTVC^2 / (s^2 + 2 zTVC wTVC s + wTVC^2)
%
%   cascaded with a Pade approximation of the 20 ms transport delay
%   (Table 1, Eq. 3). The result maps the commanded deflection delta_cmd
%   (controller output) to the actual TVC deflection delta.
%
%   Wtvc = BUILD_TVC(p, n) uses a Pade approximation of order n for the
%   delay (default n = 3). A higher order captures the destabilising phase
%   lag of the delay more faithfully near the bending frequency.
%
%   See also BUILD_NOTCH_FILTER, ASSEMBLE_LOOP.

if nargin < 2 || isempty(padeOrder), padeOrder = 3; end

s = tf('s');
Wact = p.wTVC^2 / (s^2 + 2*p.zTVC*p.wTVC*s + p.wTVC^2);

[nd, dd] = pade(p.tau, padeOrder);
Wdelay = tf(nd, dd);

Wtvc = Wact * Wdelay;
Wtvc.Name = 'TVC';
end
