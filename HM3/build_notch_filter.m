function Hx = build_notch_filter(wx, zN, zD, numSign)
% Lead-lag / notch section for bending stabilisation (Eq. 4).
%   Hx(s) = (s^2 + sgn*2*zN*wx*s + wx^2) / (s^2 + 2*zD*wx*s + wx^2)
%   INPUT
%     wx      - centre frequency [rad/s] (~ wBM +/- 4)
%     zN      - numerator damping (guideline 0.1-0.3)
%     zD      - denominator damping (guideline 0.4-0.6)
%     numSign - -1 (default, Eq. 4 as printed): RHP zeros, non-min-phase
%               phase-shaper; +1: min-phase symmetric notch
%   OUTPUT
%     Hx - tf

arguments
    wx (1,1) {mustBeNumeric, mustBeReal, mustBePositive}
    zN (1,1) {mustBeNumeric, mustBeReal, mustBeNonnegative}
    zD (1,1) {mustBeNumeric, mustBeReal, mustBePositive}
    numSign (1,1) {mustBeMember(numSign, [-1, 1])} = -1
end

num = [1, numSign*2*zN*wx, wx^2];
den = [1,         2*zD*wx, wx^2];

Hx = tf(num, den);
Hx.Name = 'Notch_Hx';
end
