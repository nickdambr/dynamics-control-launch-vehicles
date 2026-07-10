function Wtvc = build_tvc(p, padeOrder)
% TVC actuator (2nd-order) + Pade transport delay (Eq. 3): delta_cmd -> delta.
%   INPUT
%     p         - param struct (wTVC, zTVC, tau)
%     padeOrder - Pade order for the delay (default 3); higher = better phase
%                 lag near wBM
%   OUTPUT
%     Wtvc - tf, wTVC^2/(s^2+2 zTVC wTVC s+wTVC^2) * pade(tau)

arguments
    p (1,1) struct
    padeOrder (1,1) {mustBeInteger, mustBeReal, mustBePositive} = 3
end

s = tf('s');
Wact = p.wTVC^2 / (s^2 + 2*p.zTVC*p.wTVC*s + p.wTVC^2);

[nd, dd] = pade(p.tau, padeOrder);
Wdelay = tf(nd, dd);

Wtvc = Wact * Wdelay;
Wtvc.Name = 'TVC';
end
