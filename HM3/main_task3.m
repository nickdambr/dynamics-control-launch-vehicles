%% HM3 - Task 3: Robustness to parametric uncertainty (optional)
%  The aerodynamic moment coefficient mu_alpha = A6 and the control
%  effectiveness mu_c = K1 are treated as uncertain. The controller is held
%  FIXED at the Task 2 design (no re-tuning) and its stability/performance
%  robustness is assessed over the four corner cases obtained by varying
%  mu_alpha and mu_c independently by +/-30 % from nominal.
%
%  Analysis: Nichols overlay (frequency domain) + wind-gust simulations
%  (time domain) for every corner, with a summary GM/PM table.
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

%% Corner cases: [name, mu_alpha_scale, mu_c_scale]
corners = {
    'Nominal', 1.00, 1.00
    'C1',      0.70, 1.00
    'C2',      1.30, 1.00
    'C3',      1.00, 0.70
    'C4',      1.00, 1.30 };
nC = size(corners,1);

w = load_wind_profile(p0);            % same gust for all corners
L = cell(nC,1); res = cell(nC,1);
tab = zeros(nC,5);                     % [rigidGM, minGM, PM, peakTheta, peakZ]

fprintf('\n%-8s %6s %6s | %8s %8s %7s | %9s %7s %6s\n', ...
        'Case','mu_a','mu_c','rigidGM','minGM','PM','peakTh','peakZ','stab');
for i = 1:nC
    p  = load_hw3_params('mu_alpha_scale',corners{i,2},'mu_c_scale',corners{i,3});
    Gf = build_plant_full(p,'ins');
    Wf = build_tvc(p,3) * build_notch_filter(notch.wx,notch.zN,notch.zD,notch.sgn);
    [L{i}, T] = assemble_loop(Gf, K, Wf);

    am = allmargin(L{i});
    gf = am.GMFrequency; gm = 20*log10(am.GainMargin);
    rigidGM = abs(gm(find(gf>0.2 & gf<1,1)));
    minGM   = min(abs(gm));
    [~,Pm]  = margin(L{i});
    r = simulate_gust_response(T, w);
    res{i} = r;
    tab(i,:) = [rigidGM, minGM, abs(Pm), r.peak_theta*180/pi, r.peak_z];

    fprintf('%-8s %6.2f %6.2f | %7.2f  %7.2f %6.1f | %8.3f %7.2f %6d\n', ...
            corners{i,1}, corners{i,2}, corners{i,3}, ...
            rigidGM, minGM, abs(Pm), r.peak_theta*180/pi, r.peak_z, isstable(T));
end

%% ---------------------------------------------------------------- Figures
cols = lines(nC);
% Nichols overlay (phase wrapping + axis cropping for legibility)
nopt = nicholsoptions;
nopt.PhaseWrapping = 'on';
nopt.Grid = 'on';
nopt.Title.String = 'Task 3 - Nichols overlay over the +/-30% corner cases';
f1 = figure('Name','nichols_corners','Color','w','Position',[100 100 660 560]);
hN = nicholsplot(L{:}, nopt);
setoptions(hN,'XLim',[-360 0],'YLim',[-40 60]);
legend(corners(:,1),'Location','northwest');

% Gust response overlay: pitch and lateral drift
f2 = figure('Name','gust_corners','Color','w','Position',[100 100 820 360]);
tl = tiledlayout(f2,1,2,'TileSpacing','compact','Padding','compact');
title(tl,'Task 3 - Wind-gust response across the corner cases');
nexttile; hold on; grid on;
for i=1:nC, plot(res{i}.t, res{i}.theta*180/pi,'LineWidth',1.3,'Color',cols(i,:)); end
xlabel('t [s]'); ylabel('\theta [deg]'); title('Pitch attitude'); legend(corners(:,1),'Location','best');
nexttile; hold on; grid on;
for i=1:nC, plot(res{i}.t, res{i}.z,'LineWidth',1.3,'Color',cols(i,:)); end
xlabel('t [s]'); ylabel('z [m]'); title('Lateral drift'); legend(corners(:,1),'Location','best');

%% Export figures
fig_dir = fullfile(fileparts(mfilename('fullpath')), 'figures');
if ~exist(fig_dir,'dir'); mkdir(fig_dir); end
for f = [f1 f2]
    try, theme(f,'light'); catch, end
    exportgraphics(f, fullfile(fig_dir, ['task3_' get(f,'Name') '.png']), 'Resolution', 200);
end
fprintf('\nFigures written to %s\n', fig_dir);
