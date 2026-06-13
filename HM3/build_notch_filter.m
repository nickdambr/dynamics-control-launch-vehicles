function Hx = build_notch_filter(wx, zN, zD, numSign)
%BUILD_NOTCH_FILTER  Lead-lag / notch filter for bending stabilisation (Eq. 4).
%
%   Hx = BUILD_NOTCH_FILTER(wx, zN, zD) returns the second-order filter of
%   Eq. (4):
%
%       Hx(s) = (s^2 + sgn*2*zN*wx*s + wx^2) / (s^2 + 2*zD*wx*s + wx^2)
%
%   With the default sign sgn = -1 (as written in the assignment) the two
%   numerator zeros sit in the right half-plane, so the section behaves as
%   a non-minimum-phase phase-shaper: it introduces a targeted left/right
%   shift of the bending resonance on the Nichols chart, which is what
%   stabilises the lightly damped first bending mode.
%
%   Recommended ranges (assignment guidelines):
%       zN in [0.1, 0.3],  zD in [0.4, 0.6],  wx ~ wBM +/- 4 rad/s.
%
%   Hx = BUILD_NOTCH_FILTER(wx, zN, zD, +1) uses the minimum-phase
%   (symmetric notch) variant instead.
%
%   See also BUILD_TVC, ASSEMBLE_LOOP.

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
