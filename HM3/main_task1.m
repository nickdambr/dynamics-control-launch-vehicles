%% HM3 - Task 1: Rigid LV attitude control at max-qbar
%  Greensite pitch-plane model at t = 72 s (max-qbar), bending neglected,
%  ideal TVC. PD attitude controller + weak negative lateral-drift feedback,
%  tuned in frequency to |GM| ~ 6 dB / |PM| ~ 30 deg, checked on Nichols and
%  validated against a wind-gust time response.
%
%  Ref: Homework 3 (Zavoli, v1.2, May 2026), Task 1.
%  Toolboxes: Control System Toolbox (tuner uses base-MATLAB fminsearch).

clear; close all; clc;

%% Model and parameters (Table 1 @ t = 72 s)
p = load_hw3_params();
fprintf('Parameters source: %s\n', p.src);
fprintf('  A6=%.4f  K1=%.4f  V=%.1f m/s  (unstable airframe pole at +%.3f rad/s)\n', ...
        p.A6, p.K1, p.V, sqrt(p.A6));

G = build_plant_rigid(p);

%% Controller design (PD pitch + weak negative drift feedback)
fprintf('\n=== Controller tuning (rigid, ideal actuator) ===\n');
[K, m] = design_controller(G, []);    % [] => ideal actuator

% Pole-placement cross-check: the decoupled pitch CL is s^2 + K1*Kd*s + (K1*Kp - A6).
wc_eq = sqrt(p.K1*K.Kp_th - p.A6);
ze_eq = p.K1*K.Kd_th/(2*wc_eq);
fprintf(['  Equivalent pitch CL pair: wc = %.2f rad/s (course band 1-4), ' ...
         'zeta = %.2f\n'], wc_eq, ze_eq);

L = m.L;   % full open loop (drift + ideal actuator), margins classified on it
T = m.T;   % full 4-state closed loop for the gust sim

%% Rigid-body stability margins (full loop, classified by band; D'Antuono Fig 3.2)
fprintf('\n--- Rigid-body margins (full loop, classified by frequency band) ---\n');
fprintf('  Aero gain margin  : |GM| = %.2f dB @ %.2f rad/s (low-freq gain-reduction)\n', ...
        abs(m.aeroGM_dB), m.aeroGM_w);
fprintf('  Rigid phase margin: PM   = %.1f deg @ %.2f rad/s (rigid-body crossover)\n', ...
        m.rigidPM_deg, m.rigidPM_w);
fprintf('  Delay margin      : DM   = %.0f ms  (typical LV requirement >= 100 ms)\n', 1e3*m.DM_s);
fprintf('  Rigid GM / Flex margins: none (ideal actuator, no bending -> Task 2)\n');
fprintf('  Full 4-state closed loop stable (isstable): %d\n', m.stable);

%% Wind-gust time response (theta, z, zdot, delta)
w = load_wind_profile(p, Tend=80);   % 80 s horizon: dominant CL mode is tau ~ 18 s (wn=0.24, zeta=0.23), so ~5*tau to see attitude settle back to 0
r = simulate_gust_response(T, w);
fprintf('\n--- Gust response (%s gust, Vg = %.2f m/s -> peak alpha_w = %.2f deg) ---\n', ...
        w.severity, w.Vg, w.Vg/p.V*180/pi);
fprintf('  peak |theta| = %.3f deg\n', r.peak_theta*180/pi);
fprintf('  peak |z|     = %.2f m\n',   r.peak_z);
fprintf('  peak |delta| = %.3f deg\n', r.peak_delta*180/pi);
fprintf('  peak |alpha| = %.3f deg  -> peak qbar*alpha = %.1f kPa deg (qbar = %.1f kPa)\n', ...
        r.peak_alpha*180/pi, p.qbar/1000*r.peak_alpha*180/pi, p.qbar/1000);
% Lateral drift is the load-relief channel (drift-minimum), not a Nichols margin:
fprintf('  lateral drift (load-relief): peak |z| = %.1f m (<500 m), peak |zdot| = %.2f m/s (<15 m/s)\n', ...
        r.peak_z, max(abs(r.zdot)));

%% ---------------------------------------------------------------- Figures
% Full-loop Nichols in the D'Antuono Fig. 3.2 convention: the loop comes from the
% top (lateral-drift integrator), the rigid critical point sits at +180 deg, and
% the classified Aero GM / Rigid PM are marked at their crossover frequencies.
f1 = figure('Name','nichols','Color','w','Position',[100 100 680 580]);
plot_nichols_lv(L, m, 'wrange', [1e-2 1e2]);
xlim([-360 360]);
title(sprintf('Task 1 - Full-loop Nichols  (Aero |GM|=%.1f dB, Rigid PM=%.0f^\\circ)', ...
      abs(m.aeroGM_dB), m.rigidPM_deg));

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

% Angle-of-attack budget and aero load: alpha = theta + zdot/V + alpha_w
% drives qbar*alpha, the sizing quantity at max-qbar
f3 = figure('Name','alpha_load','Color','w','Position',[100 100 820 340]);
tl2 = tiledlayout(f3,1,2,'TileSpacing','compact','Padding','compact');
title(tl2,'Task 1 - Angle-of-attack budget and aerodynamic load');
nexttile; hold on; grid on;
plot(r.t, r.alpha *180/pi, 'LineWidth',1.6);
plot(r.t, r.theta *180/pi, '--','LineWidth',1.2);
plot(r.t, r.zdot/p.V*180/pi, '-.','LineWidth',1.2);
plot(r.t, r.alphaw*180/pi, ':','LineWidth',1.2);
xlabel('t [s]'); ylabel('[deg]');
legend({'$\alpha$','$\theta$','$\dot z/V$','$\alpha_w$'},'Interpreter','latex','Location','best');
title('$\alpha = \theta + \dot z/V + \alpha_w$','Interpreter','latex');
nexttile; plot(r.t, p.qbar/1000*r.alpha*180/pi, 'LineWidth',1.6); grid on;
xlabel('t [s]'); ylabel('$\bar q\,\alpha$ [kPa deg]','Interpreter','latex');
title(sprintf('Aerodynamic load  (peak %.1f kPa deg)', p.qbar/1000*r.peak_alpha*180/pi));

%% Export figures (PNG, 200 dpi) next to the script
fig_dir = fullfile(fileparts(mfilename('fullpath')), 'figures');
if ~exist(fig_dir,'dir'); mkdir(fig_dir); end
for f = [f1 f2 f3]
    try
        theme(f, 'light');    % force light theme (ignore desktop dark mode)
    catch
        f.Color = 'w';        % fallback for pre-R2025a MATLAB
    end
    exportgraphics(f, fullfile(fig_dir, ['task1_' get(f,'Name') '.png']), ...
                   'Resolution', 200);
end
fprintf('\nFigures written to %s\n', fig_dir);
