%% HM3 - Task 1: Rigid LV attitude control at max-qbar
%  Pitch-plane short-period model of the Greensite fictitious launch
%  vehicle at t = 72 s (max dynamic pressure). The bending mode is
%  neglected (rigid-body assumption) and an ideal TVC actuator is used.
%
%  A proportional-derivative attitude controller (with a weak negative
%  drift feedback on the lateral channel) is tuned in the frequency domain
%  to the assignment targets |GM| ~ 6 dB and |PM| ~ 30 deg, verified on the
%  Nichols chart, and validated through the time response to a wind gust.
%
%  Reference: Homework 3 - Attitude Control of a Launch Vehicle in
%  Atmospheric Flight (Zavoli, v1.2, May 2026), Task 1.
%  Toolboxes: Control System Toolbox (+ Optimization for the auto-tuner).

clear; close all; clc;

%% Model and parameters (Table 1 @ t = 72 s)
p = load_hw3_params();
fprintf('Parameters source: %s\n', p.src);
fprintf('  A6=%.4f  K1=%.4f  V=%.1f m/s  (unstable airframe pole at +%.3f rad/s)\n', ...
        p.A6, p.K1, p.V, sqrt(p.A6));

G = build_plant_rigid(p);

%% Controller design (PD pitch + weak negative drift feedback)
fprintf('\n=== Controller tuning (rigid, ideal actuator) ===\n');
[K, m] = design_controller(G, []);    % [] => ideal actuator Wact = 1

[L, T] = assemble_loop(G, K);

%% Frequency-domain margins
fprintf('\n--- Stability margins (Nichols / margin) ---\n');
fprintf('  Gain margin : %5.2f dB  (|GM| = %.2f dB) at w = %.3f rad/s\n', ...
        m.GM_dB, abs(m.GM_dB), m.wc_gain);
fprintf('  Phase margin: %5.1f deg (|PM| = %.1f deg) at w = %.3f rad/s\n', ...
        m.PM_deg, abs(m.PM_deg), m.wc_phase);
fprintf('  Closed-loop stable: %d\n', m.stable);

%% Wind-gust time response (theta, z, zdot, delta)
w = load_wind_profile(p);
r = simulate_gust_response(T, w);
fprintf('\n--- Gust response (%s gust, Vg = %.2f m/s) ---\n', w.severity, w.Vg);
fprintf('  peak |theta| = %.3f deg\n', r.peak_theta*180/pi);
fprintf('  peak |z|     = %.2f m\n',   r.peak_z);
fprintf('  peak |delta| = %.3f deg\n', r.peak_delta*180/pi);

%% ---------------------------------------------------------------- Figures
% Nichols chart of the open-loop transfer
f1 = figure('Name','nichols','Color','w','Position',[100 100 620 560]);
nichols(L); hold on; grid on;
ngrid;
title(sprintf('Task 1 - Rigid loop Nichols  (|GM|=%.1f dB, |PM|=%.0f^\\circ)', ...
      abs(m.GM_dB), abs(m.PM_deg)));

% Gust response: theta, z, zdot, delta
f2 = figure('Name','gust_response','Color','w','Position',[100 100 760 620]);
tl = tiledlayout(f2,2,2,'TileSpacing','compact','Padding','compact');
title(tl, sprintf('Task 1 - Rigid LV response to a %s wind gust', w.severity));

nexttile; plot(r.t, r.theta*180/pi,'LineWidth',1.4); grid on;
xlabel('t [s]'); ylabel('\theta [deg]'); title('Pitch attitude');

nexttile; plot(r.t, r.z,'LineWidth',1.4); grid on;
xlabel('t [s]'); ylabel('z [m]'); title('Lateral drift');

nexttile; plot(r.t, r.zdot,'LineWidth',1.4); grid on;
xlabel('t [s]'); ylabel('$\dot z$ [m/s]','Interpreter','latex'); title('Lateral drift rate');

nexttile; plot(r.t, r.delta*180/pi,'LineWidth',1.4); grid on;
xlabel('t [s]'); ylabel('\delta [deg]'); title('TVC deflection');

%% Export figures (PNG, 200 dpi) next to the script
fig_dir = fullfile(fileparts(mfilename('fullpath')), 'figures');
if ~exist(fig_dir,'dir'); mkdir(fig_dir); end
for f = [f1 f2]
    try, theme(f,'light'); catch, end          % publication (light) theme
    exportgraphics(f, fullfile(fig_dir, ['task1_' get(f,'Name') '.png']), ...
                   'Resolution', 200);
end
fprintf('\nFigures written to %s\n', fig_dir);
