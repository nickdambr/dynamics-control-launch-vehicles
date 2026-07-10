%% HM3 - Task 2: Full LV model (TVC dynamics, transport delay, bending mode)
%  Extends Task 1 to the full 6-state model (Eq. 1): lightly damped first
%  bending mode (wBM = 18.9 rad/s, zBM = 0.005), 2nd-order TVC (Eq. 3) + 20 ms
%  delay (Pade), and INS coupling (Eq. 2) of bending into the feedback.
%
%    Step A : rigid + PD            (warm start from Task 1)
%    Step B : + TVC + delay         (low-freq margins survive; the +39 dB
%                                     bending resonance destabilises the loop)
%    Step C : bending filter trade  (Eq.-4 lead-lag alone, deep notch, notch
%                                     triplet, notch + lead-lag)
%    Step D : wBM sensitivity       (filters fixed, true wBM perturbed +/-10 %)
%
%  Ref: Homework 3 (Zavoli, v1.2, May 2026), Task 2.
%  Toolboxes: Control System Toolbox (tuner uses base-MATLAB fminsearch).

clear; close all; clc;
warning('off','Control:analysis:MarginUnstable');

%% Model and parameters
p = load_hw3_params();

%% Step A - rigid baseline controller (Task 1 design, reused unchanged)
Grigid = build_plant_rigid(p);
fprintf('=== Step A: rigid PD design (reused from Task 1) ===\n');
[K, mR] = design_controller(Grigid, []);
[~, Trigid] = assemble_loop(Grigid, K);

%% Step B - full plant + TVC + delay, NO bending filter
Gfull = build_plant_full(p, 'ins');
Wtvc  = build_tvc(p, 3);                 % 2nd-order TVC + Pade(20 ms, order 3)
[Lb, Tb] = assemble_loop(Gfull, K, Wtvc);
fprintf('\n=== Step B: full + TVC + delay (no bending filter) ===\n');
fprintf('  |L(omega_BM)| = %.1f dB   -> closed-loop stable: %d (max Re pole = %.3g)\n', ...
        20*log10(bode(Lb,p.wBM)), isstable(Tb), max(real(pole(Tb))));

%% Step C - bending filter trade study
%  Four candidates, PD gains kept fixed (retuning usually unnecessary):
%   C-LL : Eq.-4 filter as printed (NMP numerator), swept over zN 0.1-0.3,
%          zD 0.4-0.6, wx = wBM +/- 4 rad/s. Alone it never stabilises;
%          least-unstable combo kept for the Nichols overlay.
%   C-N  : deep min-phase notch at wBM (gain stabilisation).
%   C-T  : same notch at {0.9, 1.0, 1.1} wBM - course recipe for wBM spread.
%   C-NLL: deep notch + Eq.-4 lead-lag, partner swept over the same ranges,
%          picked for max delay margin.
fprintf('\n=== Step C: bending filter trade ===\n');

wx_grid = p.wBM + (-4:2:4);
zN_grid = 0.10:0.05:0.30;
zD_grid = 0.40:0.10:0.60;

% --- C-LL: Eq.-4 lead-lag alone (least-unstable over the sweep) ---
bestLL = struct('mre',inf,'wx',NaN,'zN',NaN,'zD',NaN);
nStableLL = 0;  nTot = 0;
for wx = wx_grid
    for zN = zN_grid
        for zD = zD_grid
            nTot = nTot + 1;
            Hc = build_notch_filter(wx, zN, zD, -1);
            [~, Tc_] = assemble_loop(Gfull, K, Wtvc*Hc);
            nStableLL = nStableLL + isstable(Tc_);
            mre = max(real(pole(Tc_)));
            if mre < bestLL.mre
                bestLL = struct('mre',mre,'wx',wx,'zN',zN,'zD',zD);
            end
        end
    end
end
fprintf(['  C-LL alone: %d/%d guideline candidates stabilise the loop; ' ...
         'least unstable\n              (wx=%.1f, zN=%.2f, zD=%.2f) has max Re(pole) = %+.2f\n'], ...
        nStableLL, nTot, bestLL.wx, bestLL.zN, bestLL.zD, bestLL.mre);
Hll = build_notch_filter(bestLL.wx, bestLL.zN, bestLL.zD, -1);

% --- C-N: deep minimum-phase notch (the retained design) ---
notch.wx = p.wBM;        % rad/s (= wBM)
notch.zN = 0.002;        % numerator damping (deep, narrow null)
notch.zD = 0.7;          % denominator damping
notch.sgn = +1;          % +1 -> min-phase notch (gain stabilisation)
Hn = build_notch_filter(notch.wx, notch.zN, notch.zD, notch.sgn);

% --- C-T: course-style triplet at 0.9/1.0/1.1 wBM ---
Ht = build_notch_filter(0.9*p.wBM, notch.zN, notch.zD, +1) * ...
     build_notch_filter(    p.wBM, notch.zN, notch.zD, +1) * ...
     build_notch_filter(1.1*p.wBM, notch.zN, notch.zD, +1);

% --- C-NLL: notch + lead-lag, partner picked for max delay margin ---
bestC = struct('dm',-inf,'wx',NaN,'zN',NaN,'zD',NaN);
nStableC = 0;
for wx = wx_grid
    for zN = zN_grid
        for zD = zD_grid
            Hc = build_notch_filter(wx, zN, zD, -1);
            [Lc_, Tc_] = assemble_loop(Gfull, K, Wtvc*Hn*Hc);
            if ~isstable(Tc_), continue; end
            nStableC = nStableC + 1;
            amc = allmargin(Lc_);
            dm  = min(amc.DelayMargin);
            if dm > bestC.dm
                bestC = struct('dm',dm,'wx',wx,'zN',zN,'zD',zD);
            end
        end
    end
end
fprintf(['  C-NLL combo: %d/%d lead-lag partners stabilise the loop; ' ...
         'best (max DM):\n              wx=%.1f, zN=%.2f, zD=%.2f\n'], ...
        nStableC, nTot, bestC.wx, bestC.zN, bestC.zD);
Hcmb = Hn * build_notch_filter(bestC.wx, bestC.zN, bestC.zD, -1);

% --- comparison table over the candidates ---
cand = {'B     no filter    ', tf(1); ...
        'C-LL  lead-lag only', Hll;   ...
        'C-N   deep notch   ', Hn;    ...
        'C-T   notch triplet', Ht;    ...
        'C-NLL notch+leadlag', Hcmb};
nCand = size(cand,1);
Lcand = cell(nCand,1);
% Margins classified by frequency band (classify_margins): Aero GM (low freq),
% Rigid PM / Rigid GM (rigid body + actuator lag), and the bending attenuation
% |L(wBM)| (the deep notch gain-stabilises the mode, so there is no Flex crossover).
bands = {'w_drift',0.3*sqrt(p.A6),'w_flex',0.6*p.wBM,'w_flex_hi',1.5*p.wBM,'w_bending',p.wBM};
fprintf('\n  %-19s | %7s %7s %8s %9s %7s | %s\n', ...
        'candidate','AeroGM','RigidPM','RigidGM','|L(wBM)|','DM[ms]','stable');
for i = 1:nCand
    [Lc_, Tc_] = assemble_loop(Gfull, K, Wtvc*cand{i,2});
    Lcand{i} = minreal(Lc_, 1e-6);
    mc = classify_margins(Lcand{i}, bands{:});
    fprintf('  %-19s | %6.2f  %6.1f  %7.2f  %8.1f  %6.0f | %d\n', ...
            cand{i,1}, abs(mc.aeroGM_dB), mc.rigidPM_deg, abs(mc.rigidGM_dB), ...
            mc.LwBM_dB, 1e3*mc.DM_s, isstable(Tc_));
end

%% Step C decision - deep notch, then RE-TUNE the PD for the actuator
%  The deep notch gain-stabilises the resonance (|L(wBM)| << 0 dB) and keeps the
%  low-frequency aerodynamic gain margin. But with the Task-1 PD gains --- which
%  were designed on the IDEAL actuator --- the actuator + transport delay + notch
%  phase lag has eaten most of the rigid phase margin: it collapses from 30 deg
%  to ~15 deg, with the delay margin down at the ~100 ms guideline. The design
%  must therefore be RE-TUNED on the full loop (D'Antuono: the TVC dynamics must
%  be taken into account when placing the PD): raising the derivative gain
%  restores the phase lost to the lag, with the bending still gain-stabilised.
Wfull = Wtvc * Hn;

% --- BEFORE re-tuning: Task-1 gains (ideal-actuator design) on the full loop ---
Ktask1 = K;
[Lb1, Tb1] = assemble_loop(Gfull, Ktask1, Wfull);  Lb1 = minreal(Lb1, 1e-6);
mB = classify_margins(Lb1, bands{:});
fprintf('\n=== Deep notch, Task-1 gains (BEFORE re-tuning) ===\n');
fprintf('  Kp=%.3f Kd=%.3f | Aero |GM|=%.2f dB  Rigid PM=%.1f deg @%.2f  DM=%.0f ms | stable: %d\n', ...
        Ktask1.Kp_th, Ktask1.Kd_th, abs(mB.aeroGM_dB), mB.rigidPM_deg, mB.rigidPM_w, ...
        1e3*mB.DM_s, isstable(Tb1));

% --- RE-TUNE the PD on the full Task-2 loop (actuator + delay + notch) ---
fprintf('\n=== Re-tuned on the full loop (target Aero GM 6 dB / Rigid PM 30 deg) ===\n');
[K, mF] = design_controller(Gfull, Wfull, 'w_flex', 0.6*p.wBM, ...
              'w_flex_hi', 1.5*p.wBM, 'w_bending', p.wBM);
[Lc, Tfull] = assemble_loop(Gfull, K, Wfull);  Lc = minreal(Lc, 1e-6);
fprintf('  |L(omega_BM)| = %.1f dB (bending gain-stabilised) -> CL stable: %d (max Re pole = %.3g)\n', ...
        mF.LwBM_dB, isstable(Tfull), max(real(pole(Tfull))));
fprintf('  Aero |GM| = %.2f dB @ %.2f rad/s | Rigid GM = %.2f dB @ %.2f rad/s\n', ...
        abs(mF.aeroGM_dB), mF.aeroGM_w, abs(mF.rigidGM_dB), mF.rigidGM_w);
fprintf('  Rigid PM  = %.1f deg @ %.2f rad/s (recovered) | delay margin = %.0f ms\n', ...
        mF.rigidPM_deg, mF.rigidPM_w, 1e3*mF.DM_s);

%% Step D - sensitivity to the bending-frequency knowledge
%  Filters FIXED (designed for wBM = 18.9 rad/s); true wBM perturbed +/-10 %.
%  The deep notch needs exact wBM; the NMP lead-lag phase-stabilises the
%  resonance and buys most of the tolerance back.
fprintf('\n=== Step D: true wBM off-nominal, filters fixed ===\n');
scales = 0.90:0.05:1.10;
fprintf('  %-19s |', 'candidate');  fprintf('  x%.2f', scales);  fprintf('\n');
for i = [3 5]                       % C-N (retained) and C-NLL (robust alternative)
    fprintf('  %-19s |', cand{i,1});
    for sc = scales
        ps = load_hw3_params();  ps.wBM = sc*ps.wBM;
        Gs = build_plant_full(ps, 'ins');
        [~, Ts] = assemble_loop(Gs, K, Wtvc*cand{i,2});
        fprintf('  %5d', isstable(Ts));
    end
    fprintf('\n');
end

%% Time response to a wind gust (full model)
w  = load_wind_profile(p, Tend=80);   % 80 s horizon to match Task 1 (full-loop slow drift mode tau ~ 20 s)
rf = simulate_gust_response(Tfull, w);
rr = simulate_gust_response(Trigid, w);
fprintf('\n--- Gust response (full model, %s gust) ---\n', w.severity);
fprintf('  peak |theta| = %.3f deg, |z| = %.2f m, |delta| = %.3f deg\n', ...
        rf.peak_theta*180/pi, rf.peak_z, rf.peak_delta*180/pi);
fprintf('  peak |alpha| = %.3f deg -> peak qbar*alpha = %.1f kPa deg\n', ...
        rf.peak_alpha*180/pi, p.qbar/1000*rf.peak_alpha*180/pi);

%% ---------------------------------------------------------------- Figures
% (1) Retained design: full-loop Nichols in the launch-vehicle convention
%     (critical point +180 deg), with the classified rigid margins marked. The
%     bending mode is gain-stabilised, so the loop stays far below 0 dB near wBM
%     (no Flex crossover) rather than threading the critical point.
f1 = figure('Name','nichols','Color','w','Position',[100 100 700 600]);
plot_nichols_lv(Lc, mF, 'wrange', [1e-2 1e2], 'xlim', [-360 360], ...
    'title', sprintf(['Task 2 - Full-loop Nichols, deep notch  ' ...
             '(Aero |GM|=%.1f dB, Rigid PM=%.0f^\\circ, |L(\\omega_{BM})|=%.0f dB)'], ...
             abs(mF.aeroGM_dB), mF.rigidPM_deg, mF.LwBM_dB));

% (2) Before/after re-tuning: with the Task-1 (ideal-actuator) gains the
%     actuator + delay + notch lag collapses the rigid phase margin; re-tuning
%     the PD on the full loop restores it (same notch, +180 deg convention).
f5 = figure('Name','retune','Color','w','Position',[100 100 700 600]);
Lb1.InputName = '';  Lb1.OutputName = '';  Lc.InputName = '';  Lc.OutputName = '';
hR = nicholsplot(Lb1, Lc, {1e-2, 1e2});
setoptions(hR, 'PhaseMatching','on', 'PhaseMatchingFreq', mF.rigidPM_w, 'PhaseMatchingValue', 180, ...
              'Grid','on', 'XLim',{[-360 360]}, 'YLim',{[-40 40]}, ...
              'XLimMode','manual', 'YLimMode','manual');
title(gca, 'Task 2 - PD re-tuning on the full loop: rigid phase margin recovered');
legend({sprintf('Task-1 gains (ideal actuator): Rigid PM = %.0f^\\circ', mB.rigidPM_deg), ...
        sprintf('re-tuned on full loop: Rigid PM = %.0f^\\circ', mF.rigidPM_deg)}, ...
        'Location','southwest');

% (3) Bending-filter trade: no filter vs least-unstable lead-lag vs deep notch,
% in the launch-vehicle convention (critical point +180 deg) --- consistent with
% the other Nichols charts. A single common phase shift (from the retained deep-
% notch rigid crossover) is applied to all three loops so their bending lobes are
% directly comparable: the rigid region sits near +180 deg, the bending lobes near
% the wrapped -180 deg critical point.
f4 = figure('Name','nichols_trade','Color','w','Position',[100 100 700 600]);
ax = gca;  ngrid;  hold(ax,'on');
wv  = logspace(-2, log10(300), 5000);
[~, ph0] = bode(Lcand{3}, wv);  ph0 = squeeze(ph0);          % deep-notch reference phase
sh0 = 360*round((180 - interp1(wv, ph0, mB.rigidPM_w))/360); % common +180 shift
trcol = [0.85 0.15 0.15; 0.10 0.60 0.10; 0.10 0.20 0.85];    % red / green / blue
trnam = {'no filter (unstable)', ...
         sprintf('lead-lag Eq.4 alone (\\omega_x=%.0f, \\zeta_N=%.2f, \\zeta_D=%.1f) - marginal', ...
                 bestLL.wx, bestLL.zN, bestLL.zD), ...
         'deep notch (retained, stable)'};
ht = gobjects(3,1);
for i = 1:3
    [mag, ph] = bode(Lcand{i}, wv);  mag = squeeze(mag);  ph = squeeze(ph) + sh0;
    ht(i) = plot(ax, ph, 20*log10(mag), 'Color', trcol(i,:), 'LineWidth', 1.4, ...
                 'DisplayName', trnam{i});
end
plot(ax,  180, 0, 'r+', 'MarkerSize', 13, 'LineWidth', 1.6, 'HandleVisibility','off');
plot(ax, -180, 0, 'r+', 'MarkerSize', 13, 'LineWidth', 1.6, 'HandleVisibility','off');
xlim(ax, [-360 360]);  ylim(ax, [-40 40]);
xlabel(ax,'Open-Loop Phase (deg)');  ylabel(ax,'Open-Loop Gain (dB)');
title(ax,'Task 2 - Full-model loop: bending filter trade');
legend(ht, 'Location', 'northwest', 'FontSize', 9);

% Gust response of the full model: theta, z, zdot, delta
f2 = figure('Name','gust_response','Color','w','Position',[100 100 760 620]);
tl = tiledlayout(f2,2,2,'TileSpacing','compact','Padding','compact');
title(tl, sprintf('Task 2 - Full LV response to a %s wind gust', w.severity));
nexttile; plot(rf.t, rf.theta*180/pi,'LineWidth',1.4); grid on;
xlabel('t [s]'); ylabel('\theta [deg]'); title('Pitch attitude');
nexttile; plot(rf.t, rf.z,'LineWidth',1.4); grid on;
xlabel('t [s]'); ylabel('z [m]'); title('Lateral drift');
nexttile; plot(rf.t, rf.zdot,'LineWidth',1.4); grid on;
xlabel('t [s]'); ylabel('$\dot z$ [m/s]','Interpreter','latex'); title('Lateral drift rate');
nexttile; plot(rf.t, rf.delta*180/pi,'LineWidth',1.4); grid on;
xlabel('t [s]'); ylabel('\delta [deg]'); title('TVC deflection');

% Comparison rigid (Task 1) vs full (Task 2)
f3 = figure('Name','comparison_rigid_vs_full','Color','w','Position',[100 100 760 320]);
tl3 = tiledlayout(f3,1,3,'TileSpacing','compact','Padding','compact');
title(tl3,'Task 2 - Rigid (Task 1) vs full model response');
nexttile; plot(rr.t,rr.theta*180/pi,'--',rf.t,rf.theta*180/pi,'LineWidth',1.4); grid on;
xlabel('t [s]'); ylabel('\theta [deg]'); legend('rigid','full','Location','best');
nexttile; plot(rr.t,rr.z,'--',rf.t,rf.z,'LineWidth',1.4); grid on;
xlabel('t [s]'); ylabel('z [m]'); legend('rigid','full','Location','best');
nexttile; plot(rr.t,rr.delta*180/pi,'--',rf.t,rf.delta*180/pi,'LineWidth',1.4); grid on;
xlabel('t [s]'); ylabel('\delta [deg]'); legend('rigid','full','Location','best');

%% Export figures
fig_dir = fullfile(fileparts(mfilename('fullpath')), 'figures');
if ~exist(fig_dir,'dir'); mkdir(fig_dir); end
for f = [f1 f2 f3 f4 f5]
    try
        theme(f, 'light');    % force light theme (ignore desktop dark mode)
    catch
        f.Color = 'w';        % fallback for pre-R2025a MATLAB
    end
    exportgraphics(f, fullfile(fig_dir, ['task2_' get(f,'Name') '.png']), 'Resolution', 200);
end
fprintf('\nFigures written to %s\n', fig_dir);
