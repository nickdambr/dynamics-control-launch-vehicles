%% HM3 - Task 2: Full LV model (TVC dynamics, transport delay, bending mode)
%  The rigid design of Task 1 is extended to the full 6-state model of
%  Eq. (1): the lightly damped first bending mode (omega_BM = 18.9 rad/s,
%  zeta_BM = 0.005) is reintroduced, the ideal actuator is replaced by the
%  second-order TVC of Eq. (3) plus a 20 ms transport delay (Pade), and the
%  INS measurement model of Eq. (2) couples the bending motion into the
%  feedback. The workflow follows the assignment guidelines:
%
%    Step A : rigid + PD                (warm start from Task 1)
%    Step B : + TVC + delay             (low-freq margins survive, but the
%                                        +39 dB bending resonance destabilises
%                                        the loop)
%    Step C : bending filter trade      (four candidates: the Eq.-4 lead-lag
%                                        alone, a deep gain-stabilising notch,
%                                        a course-style notch triplet, and
%                                        notch + lead-lag)
%    Step D : wBM-knowledge sensitivity (each candidate filter held fixed
%                                        while the true bending frequency is
%                                        perturbed by up to +/-10 %)
%
%  Reference: Homework 3 (Zavoli, v1.2, May 2026), Task 2.
%  Toolboxes: Control System Toolbox (the auto-tuner uses base-MATLAB fminsearch).

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
%  Four candidates are compared before committing (PD gains kept fixed, per
%  the guideline note that retuning is usually unnecessary):
%   C-LL : the Eq.-4 filter as printed (non-minimum-phase numerator), swept
%          over the suggested ranges zN 0.1-0.3, zD 0.4-0.6, wx = wBM +/- 4
%          rad/s. Alone it never stabilises the loop: the least-unstable
%          combination is retained for the Nichols overlay.
%   C-N  : deep minimum-phase notch at wBM (gain stabilisation).
%   C-T  : triplet of the same deep notch at {0.9, 1.0, 1.1} wBM - the
%          course-notes recipe for bending-frequency uncertainty.
%   C-NLL: deep notch + Eq.-4 lead-lag ("the Notch filter and possibly other
%          filters"), partner swept over the same guideline ranges and
%          selected for max delay margin.
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
notch.wx = p.wBM;        % rad/s  (omega_BM)
notch.zN = 0.002;        % numerator damping (deep, narrow null)
notch.zD = 0.7;          % denominator damping
notch.sgn = +1;          % +1 -> minimum-phase notch (gain stabilisation)
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
fprintf('\n  %-19s | %7s %7s %7s %8s %9s | %s\n', ...
        'candidate','rigidGM','minGM','PM','DM[ms]','|L(wBM)|','stable');
for i = 1:nCand
    [Lc_, Tc_] = assemble_loop(Gfull, K, Wtvc*cand{i,2});
    Lcand{i} = Lc_;
    amc = allmargin(Lc_);
    gmc = 20*log10(amc.GainMargin);  gfc = amc.GMFrequency;
    idx = find(gfc>0.2 & gfc<1, 1);
    if isempty(idx), rgm = NaN; else, rgm = abs(gmc(idx)); end
    [~,Pmc] = margin(Lc_);
    fprintf('  %-19s | %6.2f  %6.2f  %5.1f  %7.1f  %8.1f | %d\n', ...
            cand{i,1}, rgm, min(abs(gmc)), abs(Pmc), ...
            min(amc.DelayMargin)*1000, 20*log10(bode(Lc_,p.wBM)), isstable(Tc_));
end

%% Step C decision - deep notch retained
%  The deep notch gain-stabilises the resonance (|L(wBM)| << 0 dB) while
%  staying flat at the rigid crossovers. The triplet, sized to the same
%  depth, piles up ~30 deg of phase lag at the rigid crossover (three wide
%  zD = 0.7 sections) and loses the conditionally stable low-frequency
%  margin: consistent with the course-notes philosophy, the FIRST mode is
%  too close to the rigid crossover to be blanket-notched - triplets belong
%  to the higher modes. The notch+lead-lag combo is the robust alternative
%  (see Step D) at the price of thinner nominal margins.
Wfull = Wtvc * Hn;
[Lc, Tfull] = assemble_loop(Gfull, K, Wfull);
[Gm,Pm,Wcg,Wcp] = margin(Lc);
fprintf('\n=== Retained design: deep notch (wx=%.1f, zN=%.3f, zD=%.2f) ===\n', ...
        notch.wx, notch.zN, notch.zD);
fprintf('  |L(omega_BM)| = %.1f dB   -> closed-loop stable: %d (max Re pole = %.3g)\n', ...
        20*log10(bode(Lc,p.wBM)), isstable(Tfull), max(real(pole(Tfull))));
am = allmargin(Lc);
gf = am.GMFrequency; gm = 20*log10(am.GainMargin);
idx = find(gf>0.2 & gf<1, 1);
if isempty(idx), rigidGM = NaN; else, rigidGM = gm(idx); end
fprintf('  rigid-body |GM| = %.2f dB,  |PM| = %.1f deg (preserved from Task 1)\n', ...
        abs(rigidGM), abs(Pm));
fprintf('  delay margin    = %.1f ms on top of the 20 ms already modelled\n', ...
        min(am.DelayMargin)*1000);

%% Step D - sensitivity to the bending-frequency knowledge
%  Filters stay FIXED (designed for wBM = 18.9 rad/s) while the true plant
%  bending frequency is perturbed by up to +/-10 % (the tolerance the course
%  notes cover with the notch triplet). The deep notch relies on exact
%  knowledge of wBM; adding the NMP lead-lag phase-stabilises the resonance
%  and buys most of the tolerance back.
fprintf('\n=== Step D: true wBM off-nominal, filters fixed ===\n');
scales = 0.90:0.05:1.10;
fprintf('  %-19s |', 'candidate');  fprintf('  x%.2f', scales);  fprintf('\n');
for i = 3:nCand                     % C-N, C-T, C-NLL (B and C-LL never stable)
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
w  = load_wind_profile(p);
rf = simulate_gust_response(Tfull, w);
rr = simulate_gust_response(Trigid, w);
fprintf('\n--- Gust response (full model, %s gust) ---\n', w.severity);
fprintf('  peak |theta| = %.3f deg, |z| = %.2f m, |delta| = %.3f deg\n', ...
        rf.peak_theta*180/pi, rf.peak_z, rf.peak_delta*180/pi);
fprintf('  peak |alpha| = %.3f deg -> peak qbar*alpha = %.1f kPa deg\n', ...
        rf.peak_alpha*180/pi, p.qbar/1000*rf.peak_alpha*180/pi);

%% ---------------------------------------------------------------- Figures
% Nichols: no filter vs least-unstable lead-lag vs deep notch.
% Phase wrapping + axis cropping keep the (otherwise sprawling) bending
% excursion legible around the critical region.
nopt = nicholsoptions;
nopt.PhaseWrapping = 'on';
nopt.Grid = 'on';
nopt.Title.String = 'Task 2 - Full-model loop: bending filter trade';
f1 = figure('Name','nichols','Color','w','Position',[100 100 660 560]);
hN = nicholsplot(Lcand{1},'r',Lcand{2},'g',Lcand{3},'b',nopt);
setoptions(hN,'XLim',[-360 0],'YLim',[-40 60]);
legend({'no filter (unstable)', ...
        sprintf('lead-lag Eq.4 alone (best: \\omega_x=%.0f, \\zeta_N=%.2f, \\zeta_D=%.1f) - still unstable', ...
                bestLL.wx, bestLL.zN, bestLL.zD), ...
        'deep notch (retained, stable)'}, 'Location','northwest');

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
for f = [f1 f2 f3]
    try
        theme(f, 'light');    % force light theme (ignore desktop dark mode)
    catch
        f.Color = 'w';        % fallback for pre-R2025a MATLAB
    end
    exportgraphics(f, fullfile(fig_dir, ['task2_' get(f,'Name') '.png']), 'Resolution', 200);
end
fprintf('\nFigures written to %s\n', fig_dir);
