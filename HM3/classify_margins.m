function mm = classify_margins(L, opts)
% Classify the full-loop launch-vehicle stability margins by frequency band.
%
%   Following D'Antuono (Fig. 3.2, Table 3.1) and Trotta (Table 4.1) the margins
%   of the FULL open loop L = K_FCS*G (attitude + lateral drift + actuator +
%   bending) are read at DISTINCT crossover frequencies, one per physical band:
%     Aero GM  : low-frequency aerodynamic gain-REDUCTION margin (gmdb < 0),
%                the lower edge of the conditionally-stable gain band
%     Rigid PM : phase margin at the rigid-body gain crossover (~ sqrt(A6))
%     Rigid GM : mid-frequency gain-INCREASE margin (rigid body + actuator lag;
%                absent for an ideal actuator, e.g. Task 1)
%     Flex GM/PM : gain/phase margin at the bending mode (~ w_bending)
%   Crossings below w_drift are lateral-drift artifacts (the drift position
%   integrator makes the Nichols come from the top) and are NOT reported as
%   rigid margins. Taking margin()'s default instead would pick one of these.
%
%   INPUT
%     L    - full open loop (assemble_loop), 1 + L convention
%     opts - name-value:
%              w_drift   drift/rigid boundary [rad/s] (default 0.5)
%              w_flex    rigid/flex boundary  [rad/s] (default Inf = no bending)
%              w_bending flex target freq     [rad/s] (default w_flex)
%   OUTPUT
%     mm - struct: <band>_dB / <band>_deg and matching _w frequencies
%          (NaN where a band has no crossover), DM_s, stable_am

arguments
    L {mustBeA(L, 'lti')}
    opts.w_drift   (1,1) {mustBeNumeric, mustBeReal, mustBePositive} = 0.5
    opts.w_flex    (1,1) {mustBeNumeric, mustBeReal, mustBePositive} = Inf
    opts.w_flex_hi (1,1) {mustBeNumeric, mustBeReal, mustBePositive} = Inf
    opts.w_bending (1,1) {mustBeNumeric, mustBeReal} = NaN
end
if isnan(opts.w_bending), opts.w_bending = opts.w_flex; end

warnState   = warning('off', 'Control:analysis:MarginUnstable');
restoreWarn = onCleanup(@() warning(warnState));

am   = allmargin(L);
gmdb = 20*log10(am.GainMargin(:));  gf = am.GMFrequency(:);
pm   = am.PhaseMargin(:);           pf = am.PMFrequency(:);

% Gain margins (exclude the DC / integrator entry at gf == 0)
[mm.aeroGM_dB,  mm.aeroGM_w]  = pick(gmdb, gf, gf > 0 & gf < opts.w_flex & gmdb < 0, 'minf');
[mm.rigidGM_dB, mm.rigidGM_w] = pick(gmdb, gf, gf > 0 & gf < opts.w_flex & gmdb > 0, 'minf');
[mm.flexGM_dB,  mm.flexGM_w]  = pick(gmdb, gf, gf >= opts.w_flex & gf <= opts.w_flex_hi, 'near', opts.w_bending);

% Phase margins (drift crossings below w_drift are excluded)
[mm.rigidPM_deg, mm.rigidPM_w] = pick(pm, pf, pf > opts.w_drift & pf < opts.w_flex, 'maxv');
[mm.flexPM_deg,  mm.flexPM_w]  = pick(pm, pf, pf >= opts.w_flex & pf <= opts.w_flex_hi, 'near', opts.w_bending);

% Bending-mode loop gain: for a gain-stabilised (notched) mode there is no flex
% crossover; the margin is the attenuation |L(w_bending)| below 0 dB.
if ~isnan(opts.w_bending) && isfinite(opts.w_bending)
    mm.LwBM_dB = 20*log10(abs(freqresp(L, opts.w_bending)));
else
    mm.LwBM_dB = NaN;
end

% Low-frequency lateral-drift 0 dB (gain-crossover) frequencies: conditional-
% stability crossings of the drift lobe, NOT rigid-body margins (marked as such).
mm.drift_w = pf(pf > 0 & pf <= opts.w_drift);

mm.DM_s      = min(am.DelayMargin);
mm.stable_am = am.Stable;

    function [v, w] = pick(vals, freqs, mask, mode, target)
        idx = find(mask);
        if isempty(idx), v = NaN; w = NaN; return; end
        switch mode
            case 'minf'                         % lowest crossover frequency
                [w, j] = min(freqs(idx));
            case 'near'                         % nearest a target frequency
                [~, j] = min(abs(freqs(idx) - target));  w = freqs(idx(j));
            case 'maxv'                         % largest value (rigid PM)
                [~, j] = max(vals(idx));  w = freqs(idx(j));
        end
        v = vals(idx(j));
    end
end
