function ax = plot_nichols_lv(L, mm, opts)
% Full-loop Nichols chart with the critical point at (-180 deg, 0 dB).
%
%   Uses the native MATLAB Nichols chart so the standard M/N grid with the dB
%   magnitude labels (0, +/-0.25, +/-0.5, +/-1, +/-3, +/-6, +/-12, +/-20 dB) is
%   drawn. The phase is matched so the rigid-body critical point sits at
%   -180 deg -- the course convention, from 1 + L = 0 <=> L = -1 <=>
%   (-180 deg, 0 dB) -- with the wrapped flex critical points at -540/-900.
%   (D'Antuono Fig. 3.2 shows the same chart relabeled by +360 deg; the phase
%   is defined mod 360 so the two displays are the same curve.)
%   The full open loop L = K_FCS*G comes from the top (lateral-drift position
%   integrator z/delta -> inf at DC). The classified margins (classify_margins)
%   are overlaid at their crossover frequencies, and the low-frequency drift
%   0 dB crossings are marked as conditional-stability artifacts (not margins).
%
%   INPUT
%     L    - full open loop (assemble_loop)
%     mm   - classified margins struct (classify_margins)
%     opts - name-value: wrange [wlo whi] rad/s (default [1e-3 1e3]);
%            xlim [lo hi] deg (default [-720 0]); title (default '')
%   OUTPUT
%     ax   - axis handle

arguments
    L {mustBeA(L, 'lti')}
    mm (1,1) struct
    opts.wrange (1,2) {mustBeNumeric, mustBeReal, mustBePositive} = [1e-3 1e3]
    opts.xlim   (1,2) {mustBeNumeric, mustBeReal} = [-720 0]
    opts.title  {mustBeTextScalar} = ''
end

warnState   = warning('off', 'Control:analysis:MarginUnstable');
restoreWarn = onCleanup(@() warning(warnState));

wref = mm.rigidPM_w;  if isnan(wref), wref = mm.aeroGM_w; end

% Native Nichols chart (standard M/N grid + dB labels), phase matched to -180.
L.InputName = '';  L.OutputName = '';           % drop the "from delta to delta" subtitle
h = nicholsplot(L, {opts.wrange(1), opts.wrange(2)});
setoptions(h, 'PhaseMatching','on', 'PhaseMatchingFreq', wref, 'PhaseMatchingValue', -180, ...
             'Grid','on', 'XLimMode','manual','YLimMode','manual', ...
             'XLim', {opts.xlim}, 'YLim', {[-40 40]});
ax = gca;  hold(ax, 'on');
if ~isempty(opts.title), title(ax, opts.title); end

% Overlay the classified margins, replicating the constant 360 deg phase shift.
wv = logspace(log10(opts.wrange(1)), log10(opts.wrange(2)), 4000);
[mag, ph] = bode(L, wv);  mag = squeeze(mag);  ph = squeeze(ph);  gdb = 20*log10(mag);
sh = ph + 360*round((-180 - interp1(wv, ph, wref))/360);

hleg = gobjects(0);
hleg = addmark(hleg, mm.aeroGM_w,  'Aero |GM|', 'rs');
hleg = addmark(hleg, mm.rigidPM_w, 'Rigid PM',  'rd');
hleg = addmark(hleg, mm.rigidGM_w, 'Rigid GM',  'r^');
hleg = addmark(hleg, mm.flexGM_w,  'Flex GM',   'ro');
hleg = addmark(hleg, mm.flexPM_w,  'Flex PM',   'rv');

% Lateral-drift 0 dB crossings: conditional-stability artifacts, NOT margins.
if isfield(mm, 'drift_w') && ~isempty(mm.drift_w)
    hd = gobjects(0);
    for wd = mm.drift_w(:)'
        hd = plot(ax, interp1(wv, sh, wd), interp1(wv, gdb, wd), 'kx', ...
                  'MarkerSize', 10, 'LineWidth', 1.8, 'HandleVisibility', 'off');
    end
    set(hd, 'HandleVisibility', 'on', 'DisplayName', 'drift 0 dB crossing (not a margin)');
    hleg(end+1) = hd;
end
if ~isempty(hleg)
    % Call the built-in legend on the marker handles (avoid the Nichols chart's
    % own overloaded legend method, which would treat handles as labels).
    legend(hleg, 'Location', 'southwest', 'FontSize', 8, 'Box', 'on');
end

    function hl = addmark(hl, w, name, style)
        % Mark a classified margin at frequency w (skip if absent) and add it to
        % the legend handle list with its crossover frequency.
        if isnan(w), return; end
        h2 = plot(ax, interp1(wv, sh, w), interp1(wv, gdb, w), style, ...
                  'MarkerSize', 10, 'LineWidth', 1.8, ...
                  'DisplayName', sprintf('%s (%.2g rad/s)', name, w));
        hl(end+1) = h2;
    end
end
