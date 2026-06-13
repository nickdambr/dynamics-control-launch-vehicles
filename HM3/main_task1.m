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
%  Toolboxes: Control System Toolbox (the auto-tuner uses base-MATLAB fminsearch).

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

% Cross-check against the pole-placement view of the course notes: the
% pitch-only closed loop s^2 + K1*Kd*s + (K1*Kp - A6) has
wc_eq = sqrt(p.K1*K.Kp_th - p.A6);
ze_eq = p.K1*K.Kd_th/(2*wc_eq);
fprintf(['  Equivalent pitch CL pair: wc = %.2f rad/s (course-typical 1-4), ' ...
         'zeta = %.2f (margin-driven,\n  above the 0.71-0.81 pole-placement range)\n'], wc_eq, ze_eq);

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
fprintf('\n--- Gust response (%s gust, Vg = %.2f m/s -> peak alpha_w = %.2f deg) ---\n', ...
        w.severity, w.Vg, w.Vg/p.V*180/pi);
fprintf('  peak |theta| = %.3f deg\n', r.peak_theta*180/pi);
fprintf('  peak |z|     = %.2f m\n',   r.peak_z);
fprintf('  peak |delta| = %.3f deg\n', r.peak_delta*180/pi);
fprintf('  peak |alpha| = %.3f deg  -> peak qbar*alpha = %.1f kPa deg (qbar = %.1f kPa)\n', ...
        r.peak_alpha*180/pi, p.qbar/1000*r.peak_alpha*180/pi, p.qbar/1000);

%% ---------------------------------------------------------------- Figures
% Nichols chart of the open-loop transfer
f1 = figure('Name','nichols','Color','w','Position',[100 100 620 560]);
nichols(L); hold on; grid on;
ngrid;
% Mark the points where the margins are read (phase/gain crossovers)
[magG, phG] = bode(L, m.wc_gain);    % GM: phase-crossover frequency
[magP, phP] = bode(L, m.wc_phase);   % PM: gain-crossover frequency
plot(phG(:), 20*log10(magG(:)), 'rs', 'MarkerSize',8, 'LineWidth',1.4);
plot(phP(:), 20*log10(magP(:)), 'rd', 'MarkerSize',8, 'LineWidth',1.4);
text(phG(:)+8, 20*log10(magG(:)), sprintf('|GM| @ %.2f rad/s', m.wc_gain), 'FontSize',8);
text(phP(:)+8, 20*log10(magP(:))-2, sprintf('|PM| @ %.2f rad/s', m.wc_phase), 'FontSize',8);
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

% Angle-of-attack budget and aerodynamic load indicator (course Lec. 16-17:
% alpha = theta + zdot/V + alpha_w drives the structural load qbar*alpha,
% the sizing quantity at max-qbar)
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
