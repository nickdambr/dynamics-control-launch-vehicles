%% HM3 LPV showcase — flexible vehicle: fixed vs varying notch (T008, Goal 1)
%  The frozen-time HM3 Task-2 design gain-stabilises the first bending mode
%  with a deep notch centred on omega_BM(72) = 18.9 rad/s. Over the ascent the
%  true bending frequency omega(t) sweeps 16.5 -> 31.8 rad/s, so the FIXED
%  notch detunes and eventually stops covering the resonance. This script
%  lifts the full 6-state flexible plant (BUILD_PLANT_FULL, bending + INS
%  coupling + TVC + delay) to the LPV setting and compares:
%
%    fixed   notch centred at omega(72), held over the flight  (HM3 as-is)
%    varying notch centred on omega(t)                          (LPV retune)
%
%  The LTV ode45 integration (ODE_LPV_FLEX) is the source of truth; the
%  Simulink model hm3_full_ascent_flex.slx reproduces it (RUN_FLEX_SIMULINK).
%  Reference: ticket T008. NOT part of the HM3 deliverable.

clear; close all; clc;
warning('off', 'Control:analysis:MarginUnstable');
here = fileparts(mfilename('fullpath'));
addpath(here); addpath(fileparts(here));

S  = init_simulink_lpv();
t0 = S.t0;  Tend = S.Tstop;  w72 = S.notch.wn72;

%% Frozen-time detuning sweep: |L(omega(t))| and stability, fixed vs varying
tm = (t0:5:Tend).';
[Lfix, Lvar, stbF, stbV] = deal(zeros(size(tm)));
for i = 1:numel(tm)
    p = load_hw3_params('t_ref', tm(i));
    G = build_plant_full(p, 'ins');  Wt = build_tvc(p, 3);
    [Lf, Tf] = assemble_loop(G, S.K0, Wt*build_notch_filter(w72,   0.002, 0.7, +1));
    [Lv, Tv] = assemble_loop(G, S.K0, Wt*build_notch_filter(p.wBM, 0.002, 0.7, +1));
    Lfix(i) = 20*log10(abs(squeeze(freqresp(Lf, p.wBM))));
    Lvar(i) = 20*log10(abs(squeeze(freqresp(Lv, p.wBM))));
    stbF(i) = isstable(Tf);  stbV(i) = isstable(Tv);
end
t_unstable = tm(find(~stbF, 1));            % first instant the fixed notch fails

%% Time-domain flexible response over the ascent (generator wind in the loop)
tt   = (t0:0.02:Tend).';
nt   = size(S.tvc.At, 1);
x0   = zeros(6 + 2 + nt, 1);
odeo = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
[~, xv] = ode45(@(t, x) ode_lpv_flex(t, x, make_flex(S, S.fomega)), tt, x0, odeo);  % varying
[~, xf] = ode45(@(t, x) ode_lpv_flex(t, x, make_flex(S, @(t) w72)),  tt, x0, odeo); % fixed

%% Summary
fprintf('\n=== Flexible LPV: fixed vs varying notch ===\n');
fprintf('  omega(t) sweeps %.1f -> %.1f rad/s; HM3 notch fixed at %.1f rad/s\n', ...
        min(S.fomega(tm)), max(S.fomega(tm)), w72);
fprintf('  fixed   notch: loop goes UNSTABLE at t = %g s (|L(omega)| > 0 dB)\n', t_unstable);
fprintf('  varying notch: stable all ascent, |L(omega(t))| in [%.1f, %.1f] dB\n', ...
        min(Lvar), max(Lvar));
fprintf('  peak |theta|: varying = %.3f deg (bounded) ; fixed = %.2e deg (diverges)\n', ...
        max(abs(xv(:, 3)))*180/pi, max(abs(xf(:, 3)))*180/pi);

%% ---------------------------------------------------------------- Figures
% f1 — detuning sweep: open-loop gain at the bending frequency
f1 = figure('Name', 'notch_detuning', 'Color', 'w', 'Position', [80 80 780 400]);
ax = axes(f1); hold(ax, 'on');
if ~isempty(t_unstable)
    yl = [-25 40];
    patch(ax, [t_unstable Tend Tend t_unstable], [yl(1) yl(1) yl(2) yl(2)], ...
          [0.98 0.85 0.85], 'EdgeColor', 'none', 'DisplayName', 'fixed notch unstable');
    ylim(ax, yl);
end
plot(ax, tm, Lfix, 'r.-', 'LineWidth', 1.4, 'DisplayName', 'fixed notch @ \omega(72)');
plot(ax, tm, Lvar, 'b.-', 'LineWidth', 1.4, 'DisplayName', 'varying notch @ \omega(t)');
yline(ax, 0, 'k--', '0 dB (resonance uncovered)', 'HandleVisibility', 'off');
xline(ax, 72, 'k:', 'HandleVisibility', 'off');
grid(ax, 'on'); xlabel(ax, 't [s]'); ylabel(ax, '|L(j\omega(t))|  [dB]');
title(ax, 'Open-loop gain at the bending frequency over the ascent');
legend(ax, 'Location', 'northwest');

% f2 — time-domain: attitude (clipped) and bending coordinate (log)
f2 = figure('Name', 'flex_response', 'Color', 'w', 'Position', [80 80 980 360]);
tl = tiledlayout(f2, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'Flexible-vehicle response — fixed vs varying notch');
nexttile; plot(tt, xf(:, 3)*180/pi, 'r-', tt, xv(:, 3)*180/pi, 'b-', 'LineWidth', 1.2);
ylim([-2 2]); xline(72, 'k:'); grid on; xlabel('t [s]'); ylabel('\theta [deg]');
title('Pitch attitude (\theta off-scale = divergence)');
legend('fixed', 'varying', 'Location', 'southwest');
nexttile; semilogy(tt, abs(xv(:, 5)) + 1e-12, 'b-', tt, abs(xf(:, 5)) + 1e-12, 'r-', 'LineWidth', 1.2);
xline(72, 'k:'); grid on; xlabel('t [s]'); ylabel('|\eta|  (bending coord.)');
title('Bending mode amplitude'); legend('varying', 'fixed', 'Location', 'northwest');

%% Export
fig_dir = fullfile(here, 'figures');
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end
for f = [f1 f2]
    try theme(f, 'light'); catch, set(f, 'Color', 'w'); end
    exportgraphics(f, fullfile(fig_dir, ['flex_' get(f, 'Name') '.png']), 'Resolution', 200);
end
fprintf('\nFigures written to %s\n', fig_dir);

%% ------------------------------------------------------------ local helper
function M = make_flex(S, fwn)
%MAKE_FLEX  ODE_LPV_FLEX struct with the given notch-centre handle fwn(t).
M = struct('fa1', S.fa1, 'fa3', S.fa3, 'fa4', S.fa4, 'fA6', S.fA6, 'fK1', S.fK1, ...
           'fV', S.fV, 'fomega', S.fomega, 'faqk', S.faqk, 'fsig', S.fsig, 'fphi', S.fphi, ...
           'windfun', S.windfun, 'fwn', fwn, 'zN', S.notch.zN, 'zD', S.notch.zD, ...
           'zBM', S.notch.zBM, 'At', S.tvc.At, 'Bt', S.tvc.Bt, 'Ct', S.tvc.Ct, 'Dt', S.tvc.Dt, ...
           'fKp', S.fKp, 'fKd', S.fKd, 'Kp_th0', S.K0.Kp_th, 'Kd_th0', S.K0.Kd_th, ...
           'Kp_z', S.K0.Kp_z, 'Kd_z', S.K0.Kd_z, 'sched', false);
end
