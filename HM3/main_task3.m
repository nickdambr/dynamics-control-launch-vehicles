%% HM3 - Task 3: Robustness to parametric uncertainty (optional)
%  mu_alpha = A6 (aero moment) and mu_c = K1 (control effectiveness) treated
%  as uncertain. Controller FIXED at the Task-2 design (no re-tuning),
%  assessed over the four +/-30 % box vertices V1..V4 (the assignment corner
%  cases). S1..S4 are one-at-a-time sensitivities: edge midpoints, not
%  corners, so they miss the worst combo (V3: more unstable airframe AND less
%  control authority).
%
%  Analysis: Nichols overlay + wind-gust sims per case, GM/PM/DM table.
%
%  Ref: Homework 3 (Zavoli, v1.2, May 2026), Task 3.

clear; close all; clc;
warning('off','Control:analysis:MarginUnstable');

%% Fixed controller (Task 2 retained design): PD re-tuned on the full loop + notch
%  The PD gains are the Task-2 design RE-TUNED on the full loop (actuator + delay
%  + notch), not the ideal-actuator Task-1 gains; held fixed across the corners.
p0     = load_hw3_params();
notch  = struct('wx',p0.wBM,'zN',0.002,'zD',0.7,'sgn',+1);
Gfull0 = build_plant_full(p0,'ins');
Wact0  = build_tvc(p0,3) * build_notch_filter(notch.wx,notch.zN,notch.zD,notch.sgn);
K = design_controller(Gfull0, Wact0, 'w_flex',0.6*p0.wBM, 'w_flex_hi',1.5*p0.wBM, ...
                      'w_bending',p0.wBM, 'verbose',false);
fprintf('Fixed controller (Task-2 re-tuned): Kp_th=%.3f Kd_th=%.3f | notch wx=%.1f zN=%.3f zD=%.2f\n', ...
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

w = load_wind_profile(p0, Tend=80);    % same gust for all cases; 80 s horizon to match Task 1
L = cell(nC,1); res = cell(nC,1); mm = cell(nC,1);

fprintf('\n%-8s %6s %6s | %7s %7s %8s %7s | %8s %7s %5s\n', ...
        'Case','mu_a','mu_c','AeroGM','RigidPM','RigidGM','DM[ms]','peakTh','peakZ','stab');
for i = 1:nC
    p  = load_hw3_params('mu_alpha_scale',cases{i,2},'mu_c_scale',cases{i,3});
    Gf = build_plant_full(p,'ins');
    Wf = build_tvc(p,3) * build_notch_filter(notch.wx,notch.zN,notch.zD,notch.sgn);
    [L{i}, T] = assemble_loop(Gf, K, Wf);
    L{i} = minreal(L{i}, 1e-6);

    % margins classified per corner (aero-pole band follows the corner's A6)
    mm{i} = classify_margins(L{i}, 'w_drift',0.3*sqrt(p.A6), 'w_flex',0.6*p.wBM, ...
                             'w_flex_hi',1.5*p.wBM, 'w_bending',p.wBM);
    r = simulate_gust_response(T, w);
    res{i} = r;

    fprintf('%-8s %6.2f %6.2f | %6.2f  %6.1f  %7.2f  %6.0f | %7.3f %6.2f %5d\n', ...
            cases{i,1}, cases{i,2}, cases{i,3}, ...
            abs(mm{i}.aeroGM_dB), mm{i}.rigidPM_deg, abs(mm{i}.rigidGM_dB), 1e3*mm{i}.DM_s, ...
            r.peak_theta*180/pi, r.peak_z, isstable(T));
end
fprintf(['\n(V* = uncertainty-box vertices, the assignment corner cases;' ...
         ' S* = one-at-a-time sensitivities. AeroGM/RigidGM in dB, RigidPM in deg.)\n']);

%% ---------------------------------------------------------------- Figures
cols = lines(nPlot);
% Nichols overlay over Nominal + vertices, critical point at (-180 deg, 0 dB)
% (course convention, from 1 + L = 0). A single common phase shift (from the
% nominal rigid crossover) is applied to all corners so their spread is directly
% comparable. Zoomed on the rigid region: V3 (max instability, min authority) is
% the corner whose aerodynamic gain margin nearly vanishes -- its curve rides
% closest to the -180 critical point.
f1 = figure('Name','nichols_corners','Color','w','Position',[100 100 700 600]);
ax = gca;  ngrid;  hold(ax,'on');
% ngrid draws its M/N grid on the [0,360] sheet; move it onto [-360,0]
% (critical point at -180) so it backs the plotted window
for k = 1:numel(ax.Children)
    hg = ax.Children(k);
    if isprop(hg, 'XData'), hg.XData = hg.XData - 360;
    else, hg.Position(1) = hg.Position(1) - 360; end
end
wv  = logspace(-2, log10(30), 3000);
[~, ph1] = bode(L{1}, wv);  ph1 = squeeze(ph1);       % nominal phase (deg, unwrapped)
sh0 = 360*round((-180 - interp1(wv, ph1, mm{1}.rigidPM_w))/360);   % common -180 shift
hc = gobjects(nPlot,1);
for i = 1:nPlot
    [mag, ph] = bode(L{i}, wv);  mag = squeeze(mag);  ph = squeeze(ph) + sh0;
    hc(i) = plot(ax, ph, 20*log10(mag), 'Color', cols(i,:), 'LineWidth', 1.5, ...
                 'DisplayName', cases{i,1});
end
plot(ax, -180, 0, 'r+', 'MarkerSize', 13, 'LineWidth', 1.6, 'HandleVisibility','off');
xlim(ax, [-270 -90]);  ylim(ax, [-15 20]);
xlabel(ax,'Open-Loop Phase (deg)');  ylabel(ax,'Open-Loop Gain (dB)');
title(ax,'Task 3 - Nichols overlay over the \pm30% corner cases');
legend(hc, 'Location', 'southwest', 'FontSize', 9);

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
    try
        theme(f, 'light');    % force light theme (ignore desktop dark mode)
    catch
        f.Color = 'w';        % fallback for pre-R2025a MATLAB
    end
    exportgraphics(f, fullfile(fig_dir, ['task3_' get(f,'Name') '.png']), 'Resolution', 200);
end
fprintf('\nFigures written to %s\n', fig_dir);
