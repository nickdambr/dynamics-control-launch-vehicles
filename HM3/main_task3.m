%% HM3 - Task 3: Robustness to parametric uncertainty (optional)
%  The aerodynamic moment coefficient mu_alpha = A6 and the control
%  effectiveness mu_c = K1 are treated as uncertain. The controller is held
%  FIXED at the Task 2 design (no re-tuning) and its stability/performance
%  robustness is assessed over the FOUR CORNER CASES of the uncertainty box,
%  i.e. the vertices obtained by independently varying mu_alpha and mu_c by
%  +/-30 % (V1..V4). The four one-at-a-time variations (S1..S4) are kept as
%  a complementary sensitivity study: they attribute the effects to each
%  parameter, but they sit on the edge midpoints of the box, not on its
%  corners, and therefore miss the worst-case combination (V3: more unstable
%  airframe AND less control authority).
%
%  Analysis: Nichols overlay (frequency domain) + wind-gust simulations
%  (time domain) for every case, with a summary GM/PM/delay-margin table.
%
%  Reference: Homework 3 (Zavoli, v1.2, May 2026), Task 3.

clear; close all; clc;
warning('off','Control:analysis:MarginUnstable');

%% Fixed controller (Task 2 design): rigid PD + bending notch
p0     = load_hw3_params();
Grigid = build_plant_rigid(p0);
K      = design_controller(Grigid, [], 'verbose', false);
notch  = struct('wx',p0.wBM,'zN',0.002,'zD',0.7,'sgn',+1);
fprintf('Fixed controller: Kp_th=%.3f Kd_th=%.3f | notch wx=%.1f zN=%.3f zD=%.2f\n', ...
        K.Kp_th, K.Kd_th, notch.wx, notch.zN, notch.zD);

%% Cases: box vertices (corner cases) + one-at-a-time sensitivities
%  [name, mu_alpha_scale, mu_c_scale]
cases = {
    'Nominal', 1.00, 1.00
    'V1',      0.70, 0.70      % box vertices = the assignment corner cases
    'V2',      0.70, 1.30
    'V3',      1.30, 0.70      % worst case: max instability, min authority
    'V4',      1.30, 1.30
    'S1',      0.70, 1.00      % one-at-a-time sensitivities (extra)
    'S2',      1.30, 1.00
    'S3',      1.00, 0.70
    'S4',      1.00, 1.30 };
nC    = size(cases,1);
nPlot = 5;                             % figures show Nominal + V1..V4 only

w = load_wind_profile(p0);             % same gust for all cases
L = cell(nC,1); res = cell(nC,1);

fprintf('\n%-8s %6s %6s | %8s %8s %7s %8s | %9s %7s %6s\n', ...
        'Case','mu_a','mu_c','rigidGM','minGM','PM','DM[ms]','peakTh','peakZ','stab');
for i = 1:nC
    p  = load_hw3_params('mu_alpha_scale',cases{i,2},'mu_c_scale',cases{i,3});
    Gf = build_plant_full(p,'ins');
    Wf = build_tvc(p,3) * build_notch_filter(notch.wx,notch.zN,notch.zD,notch.sgn);
    [L{i}, T] = assemble_loop(Gf, K, Wf);

    am = allmargin(L{i});
    gf = am.GMFrequency; gm = 20*log10(am.GainMargin);
    rigidGM = abs(gm(find(gf>0.2 & gf<1,1)));
    minGM   = min(abs(gm));
    dm      = min(am.DelayMargin);
    [~,Pm]  = margin(L{i});
    r = simulate_gust_response(T, w);
    res{i} = r;

    fprintf('%-8s %6.2f %6.2f | %7.2f  %7.2f %6.1f %8.1f | %8.3f %7.2f %6d\n', ...
            cases{i,1}, cases{i,2}, cases{i,3}, ...
            rigidGM, minGM, abs(Pm), dm*1000, ...
            r.peak_theta*180/pi, r.peak_z, isstable(T));
end
fprintf(['\n(V* = uncertainty-box vertices, the assignment corner cases;' ...
         ' S* = one-at-a-time sensitivities.)\n']);

%% ---------------------------------------------------------------- Figures
cols = lines(nPlot);
% Nichols overlay over Nominal + vertices (phase wrapping + cropping)
nopt = nicholsoptions;
nopt.PhaseWrapping = 'on';
nopt.Grid = 'on';
nopt.Title.String = 'Task 3 - Nichols overlay over the +/-30% corner cases';
f1 = figure('Name','nichols_corners','Color','w','Position',[100 100 660 560]);
hN = nicholsplot(L{1:nPlot}, nopt);
setoptions(hN,'XLim',[-360 0],'YLim',[-40 60]);
legend(cases(1:nPlot,1),'Location','northwest');

% Gust response overlay: pitch and lateral drift (Nominal + vertices)
f2 = figure('Name','gust_corners','Color','w','Position',[100 100 820 360]);
tl = tiledlayout(f2,1,2,'TileSpacing','compact','Padding','compact');
title(tl,'Task 3 - Wind-gust response across the corner cases');
nexttile; hold on; grid on;
for i=1:nPlot, plot(res{i}.t, res{i}.theta*180/pi,'LineWidth',1.3,'Color',cols(i,:)); end
xlabel('t [s]'); ylabel('\theta [deg]'); title('Pitch attitude');
legend(cases(1:nPlot,1),'Location','best');
nexttile; hold on; grid on;
for i=1:nPlot, plot(res{i}.t, res{i}.z,'LineWidth',1.3,'Color',cols(i,:)); end
xlabel('t [s]'); ylabel('z [m]'); title('Lateral drift');
legend(cases(1:nPlot,1),'Location','best');

%% Export figures
fig_dir = fullfile(fileparts(mfilename('fullpath')), 'figures');
if ~exist(fig_dir,'dir'); mkdir(fig_dir); end
for f = [f1 f2]
    try, theme(f,'light'); catch, end
    exportgraphics(f, fullfile(fig_dir, ['task3_' get(f,'Name') '.png']), 'Resolution', 200);
end
fprintf('\nFigures written to %s\n', fig_dir);
