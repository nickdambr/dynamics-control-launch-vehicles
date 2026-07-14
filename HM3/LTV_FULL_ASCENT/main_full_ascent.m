%% HM3 full-ascent LPV attitude control (beyond the assignment)
%  Extends the frozen-time HM3 design (plant/PD/notch/TVC at max-q, t_ref=72 s)
%  to the whole ascent. The rigid pitch-plane plant (BUILD_PLANT_RIGID) is made
%  time-varying from GreensiteLPV_DATA.mat, so the wind generator
%  (strong_wind.slx, 0-140 s) and the dynamics share one clock. Compares on the
%  same LTV plant:
%    frozen    - single max-q PD held over the flight (shows where it degrades)
%    scheduled - PD gains Kp_th(t), Kd_th(t) from DESIGN_CONTROLLER on a grid
%                of frozen plants (continuation)
%
%  ode45 here is the source of truth; hm3_full_ascent.slx
%  (BUILD_HM3_FULL_ASCENT) reproduces it (RUN_FULL_ASCENT_SIMULINK). Portfolio
%  showcase, NOT part of the HM3 deliverable (assignment is the max-q point only).
%  Reference: Homework 3 (Zavoli, v1.2, May 2026); ticket T007.
%  Toolboxes: Control System Toolbox (schedule reuses HM3's fminsearch tuner).

clear; close all; clc;
warning('off', 'Control:analysis:MarginUnstable');
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));                       % HM3 helpers (load_hw3_params, ...)

%% Setup: LPV coefficients, gain schedule, full-ascent wind (pushed to base too)
S  = init_simulink_lpv();
t0 = S.t0;  Tend = S.Tstop;

%% LTV closed-loop integration: frozen gains vs gain schedule
tt    = (t0:0.02:Tend).';
x0    = zeros(4, 1);                             % [z, zdot, theta, thetadot]
odeo  = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

[~, xF] = ode45(@(t, x) ode_lpv_ascent(t, x, make_model(S, 0)), tt, x0, odeo);
[~, xS] = ode45(@(t, x) ode_lpv_ascent(t, x, make_model(S, 1)), tt, x0, odeo);

rF = unpack(tt, xF, S, 0);                       % frozen-gain response
rS = unpack(tt, xS, S, 1);                       % gain-scheduled response

%% Frozen-time margin sweep along the trajectory
%  Freeze the plant at each time, read loop margins with the fixed max-q gains
%  (how far the single-point design stretches) and with the scheduled gains.
tm = (t0:2.5:Tend).';
[gmF, pmF, gmS, pmS] = deal(zeros(size(tm)));
for i = 1:numel(tm)
    Gi = build_plant_rigid(load_hw3_params('t_ref', tm(i)));
    [gmF(i), pmF(i)] = loop_margin(Gi, S.K0);                 % frozen gains
    Ksi = S.K0; Ksi.Kp_th = S.fKp(tm(i)); Ksi.Kd_th = S.fKd(tm(i));
    [gmS(i), pmS(i)] = loop_margin(Gi, Ksi);                  % scheduled gains
end

%% Consistency check at t_ref: the LPV loop must reduce to the HM3 frozen loop
p72 = load_hw3_params();                          % t_ref = 72 s
[~, T72] = assemble_loop(build_plant_rigid(p72), S.K0, []);
wg72 = load_wind_profile(p72);                    % HM3 1-cosine gust (12 s)
rHM3 = simulate_gust_response(T72, wg72);          % HM3 frozen-time response
Mf = make_model(S, 0);                             % LPV model, coeffs frozen at 72 s
cc = @(v) griddedInterpolant([0 200], [v v], 'linear', 'nearest');
Mf.fc1 = cc(S.fc1(72)); Mf.fc2 = cc(S.fc2(72)); Mf.fc3 = cc(S.fc3(72));
Mf.fc4 = cc(S.fc4(72)); Mf.fc5 = cc(S.fc5(72)); Mf.fc6 = cc(S.fc6(72));
Mf.fc7 = cc(S.fc7(72));
Mf.windfun = griddedInterpolant(wg72.t(:), wg72.alphaw(:), 'linear', 'nearest');
[~, xf72] = ode45(@(t, x) ode_lpv_ascent(t, x, Mf), wg72.t, x0, odeo);
err_consistency = max(abs(xf72(:, 3) - rHM3.theta));

%% Summary
fprintf('\n=== Full-ascent LTV summary (rigid LPV plant, %g-%g s) ===\n', t0, Tend);
fprintf('  %-12s | peak|theta|  peak|z|  peak|delta|  peak qbar*alpha\n', 'controller');
fprintf('  %-12s |  %7.3f deg  %6.2f m  %8.3f deg  %9.1f kPa.deg\n', 'frozen', ...
        rF.pk_theta, rF.pk_z, rF.pk_delta, rF.pk_qa);
fprintf('  %-12s |  %7.3f deg  %6.2f m  %8.3f deg  %9.1f kPa.deg\n', 'scheduled', ...
        rS.pk_theta, rS.pk_z, rS.pk_delta, rS.pk_qa);
fprintf('  consistency @ t_ref=72 s: max|theta_LPV - theta_HM3| = %.2e rad\n', err_consistency);
[qmax, iqm] = max(S.Q);  tq = S.tg(iqm);
fprintf('  dataset max-q = %.1f kPa at t = %.0f s  (HM3 design point t_ref = 72 s)\n', ...
        qmax/1000, tq);

%% ---------------------------------------------------------------- Figures
% f1 — closed-loop response: frozen vs scheduled
f1 = figure('Name', 'response', 'Color', 'w', 'Position', [80 80 980 320]);
tl1 = tiledlayout(f1, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl1, 'Full-ascent LTV response — frozen vs gain-scheduled PD');
nexttile; plot(rF.t, rF.theta*180/pi, 'b-', rS.t, rS.theta*180/pi, 'r-', 'LineWidth', 1.2);
xline(72, 'k:'); grid on; xlabel('t [s]'); ylabel('\theta [deg]');
title('Pitch attitude'); legend('frozen', 'scheduled', 'Location', 'best');
nexttile; plot(rF.t, rF.z, 'b-', rS.t, rS.z, 'r-', 'LineWidth', 1.2);
xline(72, 'k:'); grid on; xlabel('t [s]'); ylabel('z [m]'); title('Lateral drift');
nexttile; plot(rF.t, rF.delta*180/pi, 'b-', rS.t, rS.delta*180/pi, 'r-', 'LineWidth', 1.2);
xline(72, 'k:'); grid on; xlabel('t [s]'); ylabel('\delta [deg]'); title('TVC command');

% f2 — structural-load indicator qbar*alpha (what the max-q design protects)
f2 = figure('Name', 'qbar_alpha', 'Color', 'w', 'Position', [80 80 720 360]);
ax2 = axes(f2); hold(ax2, 'on');
yl = max([rF.qa; rS.qa]) * 1.1;
patch(ax2, [tq-10 tq+10 tq+10 tq-10], [-yl -yl yl yl], [0.95 0.9 0.7], ...
      'EdgeColor', 'none', 'FaceAlpha', 0.6, 'DisplayName', 'max-q region');
plot(ax2, rF.t, rF.qa, 'b-', 'LineWidth', 1.3, 'DisplayName', 'frozen');
plot(ax2, rS.t, rS.qa, 'r-', 'LineWidth', 1.3, 'DisplayName', 'scheduled');
xline(ax2, 72, 'k:', 'HandleVisibility', 'off');
grid(ax2, 'on'); xlabel(ax2, 't [s]'); ylabel(ax2, '$\bar{q}\,\alpha$ [kPa$\cdot$deg]', 'Interpreter', 'latex');
title(ax2, 'Structural-load indicator over the ascent'); legend(ax2, 'Location', 'best');

% f3 — frozen-time margin sweep
f3 = figure('Name', 'margin_sweep', 'Color', 'w', 'Position', [80 80 760 540]);
tl3 = tiledlayout(f3, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl3, 'Frozen-time loop margins along the trajectory');
nexttile; plot(tm, abs(gmF), 'b.-', tm, abs(gmS), 'r.-', 'LineWidth', 1.2);
yline(6, 'k--', 'target 6 dB'); xline(72, 'k:'); grid on;
xlabel('t [s]'); ylabel('|GM| [dB]'); legend('frozen gains', 'scheduled', 'Location', 'best');
title('Gain margin (magnitude)');
nexttile; plot(tm, abs(pmF), 'b.-', tm, abs(pmS), 'r.-', 'LineWidth', 1.2);
yline(30, 'k--', 'target 30^\circ'); xline(72, 'k:'); grid on;
xlabel('t [s]'); ylabel('|PM| [deg]'); legend('frozen gains', 'scheduled', 'Location', 'best');
title('Phase margin (magnitude)');

%% Export
fig_dir = fullfile(here, 'figures');
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end
for f = [f1 f2 f3]
    try theme(f, 'light'); catch, set(f, 'Color', 'w'); end
    exportgraphics(f, fullfile(fig_dir, ['fullascent_' get(f, 'Name') '.png']), 'Resolution', 200);
end
fprintf('\nFigures written to %s\n', fig_dir);

%% ------------------------------------------------------------ local helpers
function M = make_model(S, sched)
% Pack INIT_SIMULINK_LPV data into the ODE_LPV_ASCENT struct.
%   INPUT
%     S     - init_simulink_lpv struct
%     sched - 0 frozen gains, 1 scheduled
%   OUTPUT
%     M - struct for ODE_LPV_ASCENT
M = struct('fc1', S.fc1, 'fc2', S.fc2, 'fc3', S.fc3, 'fc4', S.fc4, ...
           'fc5', S.fc5, 'fc6', S.fc6, 'fc7', S.fc7, 'windfun', S.windfun, ...
           'fKp', S.fKp, 'fKd', S.fKd, 'Kp_th0', S.K0.Kp_th, 'Kd_th0', S.K0.Kd_th, ...
           'Kp_z', S.K0.Kp_z, 'Kd_z', S.K0.Kd_z, 'sched', logical(sched));
end

function r = unpack(tt, x, S, sched)
% Time histories, load indicator and peak metrics from an ode45 solution.
% delta reconstructed with the gains used in this run (frozen or scheduled).
%   INPUT
%     tt    - time grid
%     x     - ode45 states (Nx4)
%     S     - init_simulink_lpv struct
%     sched - 0 frozen gains, 1 scheduled
%   OUTPUT
%     r - struct: t, theta, z, zdot, thetadot, delta, alpha, qa, peaks
if sched
    Kp = S.fKp(tt);  Kd = S.fKd(tt);
else
    Kp = S.K0.Kp_th*ones(size(tt));  Kd = S.K0.Kd_th*ones(size(tt));
end
delta = -(Kp.*x(:, 3) + Kd.*x(:, 4) + S.K0.Kp_z*x(:, 1) + S.K0.Kd_z*x(:, 2));
V     = S.fV(tt);
alpha = x(:, 3) + x(:, 2)./V - S.windfun(tt);    % total angle of attack
                                                 % (minus: same convention as the
                                                 % RHS, ode_lpv_ascent.m lines 26/28)
qa    = (S.fQ(tt)/1000) .* (alpha*180/pi);       % qbar*alpha [kPa.deg]
r = struct('t', tt, 'theta', x(:, 3), 'z', x(:, 1), 'zdot', x(:, 2), ...
           'thetadot', x(:, 4), 'delta', delta, 'alpha', alpha, 'qa', qa);
r.pk_theta = max(abs(r.theta))*180/pi;
r.pk_z     = max(abs(r.z));
r.pk_delta = max(abs(r.delta))*180/pi;
r.pk_qa    = max(abs(r.qa));
end

function [gm_db, pm_deg] = loop_margin(G, K)
% Gain/phase margin of the rigid loop (ideal actuator).
%   INPUT
%     G - plant
%     K - gains struct
%   OUTPUT
%     gm_db  - gain margin [dB]
%     pm_deg - phase margin [deg]
[L, ~] = assemble_loop(G, K, []);
[g, p] = margin(L);
gm_db = 20*log10(g);  pm_deg = p;
end
