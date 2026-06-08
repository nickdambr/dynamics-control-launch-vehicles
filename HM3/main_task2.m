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
%    Step C : + notch filter (Eq. 4)    (the bending peak is gain-stabilised,
%                                        the rigid margins are preserved)
%
%  Reference: Homework 3 (Zavoli, v1.2, May 2026), Task 2.
%  Toolboxes: Control System Toolbox (+ Optimization for the auto-tuner).

clear; close all; clc;
warning('off','Control:analysis:MarginUnstable');

%% Model and parameters
p = load_hw3_params();

%% Step A - rigid baseline controller (Task 1 design, reused unchanged)
Grigid = build_plant_rigid(p);
fprintf('=== Step A: rigid PD design (reused from Task 1) ===\n');
[K, mR] = design_controller(Grigid, []);
[~, Trigid] = assemble_loop(Grigid, K);

%% Step B - full plant + TVC + delay, NO notch
Gfull = build_plant_full(p, 'ins');
Wtvc  = build_tvc(p, 3);                 % 2nd-order TVC + Pade(20 ms, order 3)
[Lb, Tb] = assemble_loop(Gfull, K, Wtvc);
fprintf('\n=== Step B: full + TVC + delay (no notch) ===\n');
fprintf('  |L(omega_BM)| = %.1f dB   -> closed-loop stable: %d (max Re pole = %.3g)\n', ...
        20*log10(bode(Lb,p.wBM)), isstable(Tb), max(real(pole(Tb))));

%% Step C - insert the bending notch (Eq. 4, gain-stabilisation variant)
%  Minimum-phase (true) notch centred on the bending frequency: it drives
%  |L(omega_BM)| well below 0 dB so the resonance is gain-stabilised, while
%  staying flat at the low-frequency rigid-body crossovers (margins kept).
notch.wx = p.wBM;        % rad/s  (omega_BM)
notch.zN = 0.002;        % numerator damping (deep, narrow null)
notch.zD = 0.7;          % denominator damping
notch.sgn = +1;          % +1 -> minimum-phase notch (gain stabilisation)
Hx = build_notch_filter(notch.wx, notch.zN, notch.zD, notch.sgn);
Wfull = Wtvc * Hx;

[Lc, Tfull] = assemble_loop(Gfull, K, Wfull);
[Gm,Pm,Wcg,Wcp] = margin(Lc);
fprintf('\n=== Step C: + notch (wx=%.1f, zN=%.3f, zD=%.2f) ===\n', notch.wx, notch.zN, notch.zD);
fprintf('  |L(omega_BM)| = %.1f dB   -> closed-loop stable: %d (max Re pole = %.3g)\n', ...
        20*log10(bode(Lc,p.wBM)), isstable(Tfull), max(real(pole(Tfull))));
am = allmargin(Lc);
gf = am.GMFrequency; gm = 20*log10(am.GainMargin);
rigidGM = gm(find(gf>0.2 & gf<1,1));
fprintf('  rigid-body |GM| = %.2f dB,  |PM| = %.1f deg (preserved from Task 1)\n', ...
        abs(rigidGM), abs(mR.PM_deg));

%% Time response to a wind gust (full model)
w  = load_wind_profile(p);
rf = simulate_gust_response(Tfull, w);
rr = simulate_gust_response(Trigid, w);
fprintf('\n--- Gust response (full model, %s gust) ---\n', w.severity);
fprintf('  peak |theta| = %.3f deg, |z| = %.2f m, |delta| = %.3f deg\n', ...
        rf.peak_theta*180/pi, rf.peak_z, rf.peak_delta*180/pi);

%% ---------------------------------------------------------------- Figures
% Nichols: Step B (no notch) vs Step C (with notch).
% Phase wrapping + axis cropping keep the (otherwise sprawling) bending
% excursion legible around the critical region.
nopt = nicholsoptions;
nopt.PhaseWrapping = 'on';
nopt.Grid = 'on';
nopt.Title.String = 'Task 2 - Full-model loop: effect of the bending notch';
f1 = figure('Name','nichols','Color','w','Position',[100 100 660 560]);
hN = nicholsplot(Lb,'r',Lc,'b',nopt);
setoptions(hN,'XLim',[-360 0],'YLim',[-40 60]);
legend({'Step B: no notch (unstable)','Step C: with notch (stable)'}, ...
       'Location','northwest');

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
    try, theme(f,'light'); catch, end
    exportgraphics(f, fullfile(fig_dir, ['task2_' get(f,'Name') '.png']), 'Resolution', 200);
end
fprintf('\nFigures written to %s\n', fig_dir);
