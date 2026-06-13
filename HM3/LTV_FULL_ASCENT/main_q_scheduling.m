%% HM3 LPV showcase — scheduling on q(t) instead of t (T008, Goal 2)
%  The gain schedule of MAIN_FULL_ASCENT is keyed on flight TIME. A flight
%  controller cannot measure time-since-launch directly; the textbook LPV
%  choice is to schedule on a measurable parameter — here the dynamic
%  pressure q (the quantity that drives the aerodynamic instability). This
%  script re-keys the same DESIGN_CONTROLLER gains on q(t) and asks the honest
%  question: is q a good scheduling variable for this vehicle?
%
%  Finding: NO, not cleanly. The aerodynamic instability A6(t) peaks at
%  t ~ 72 s but the dynamic pressure Q(t) peaks earlier (~65-67 s), so the
%  gain-vs-q map is HYSTERETIC — the ascending and descending q branches need
%  different gains at the same q — and Q plateaus late in flight while the
%  gains keep falling. The script quantifies the resulting schedule error.
%  (Mach, which is monotonic, would be the better measurable; noted in README.)
%
%  Reference: ticket T008. Rigid LPV plant (Goal 2 is independent of bending).

clear; close all; clc;
warning('off', 'Control:analysis:MarginUnstable');
here = fileparts(mfilename('fullpath'));
addpath(here); addpath(fileparts(here));

S  = init_simulink_lpv();
t0 = S.t0;  Tend = S.Tstop;

%% Build a q-keyed gain lookup from the ASCENDING branch (q monotonic there)
Qs = S.fQ(S.tsched)/1000;                 % dyn. pressure at schedule points [kPa]
[~, ipk] = max(Qs);                       % q-peak splits ascending / descending
qa  = Qs(1:ipk);  Kpa = S.Kp_sched(1:ipk);  Kda = S.Kd_sched(1:ipk);
Kp_of_q = griddedInterpolant(qa, Kpa, 'linear', 'nearest');
Kd_of_q = griddedInterpolant(qa, Kda, 'linear', 'nearest');
% gains as a function of TIME when scheduled on the measured q(t)
fKp_q = @(t) Kp_of_q(S.fQ(t)/1000);
fKd_q = @(t) Kd_of_q(S.fQ(t)/1000);

% hysteresis: gain the q-lookup commands on the descending branch vs the
% gain that branch actually needs (the t-schedule), at matched q
Kp_cmd_dn = fKp_q(S.tsched(ipk+1:end));
Kp_need_dn = S.Kp_sched(ipk+1:end);
hyst = max(abs(Kp_cmd_dn - Kp_need_dn) ./ Kp_need_dn) * 100;

%% Closed-loop responses: frozen vs t-scheduled vs q-scheduled
tt   = (t0:0.02:Tend).';
odeo = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
runs = struct('name', {'frozen', 't-sched', 'q-sched'});
[runs.r] = deal([]);
mdls = {make_model(S, 'frozen', fKp_q, fKd_q), ...
        make_model(S, 't',      fKp_q, fKd_q), ...
        make_model(S, 'q',      fKp_q, fKd_q)};
for k = 1:3
    [~, x] = ode45(@(t, x) ode_lpv_ascent(t, x, mdls{k}), tt, zeros(4, 1), odeo);
    runs(k).r = unpack(tt, x, S, mdls{k});
end

%% Frozen-time margin sweep: t-scheduled vs q-scheduled gains
tm = (t0:2.5:Tend).';
[gmT, pmT, gmQ, pmQ] = deal(zeros(size(tm)));
for i = 1:numel(tm)
    Gi = build_plant_rigid(load_hw3_params('t_ref', tm(i)));
    Kt = S.K0; Kt.Kp_th = S.fKp(tm(i)); Kt.Kd_th = S.fKd(tm(i));
    [gmT(i), pmT(i)] = loop_margin(Gi, Kt);
    Kq = S.K0; Kq.Kp_th = fKp_q(tm(i)); Kq.Kd_th = fKd_q(tm(i));
    [gmQ(i), pmQ(i)] = loop_margin(Gi, Kq);
end

%% Summary
fprintf('\n=== q- vs t-scheduling (rigid LPV plant) ===\n');
fprintf('  Q peaks at t=%g s; A6 (instability) peaks at t~72 s -> hysteresis\n', S.tsched(ipk));
fprintf('  max gain mismatch of the q-lookup on the descending branch: %.0f %%\n', hyst);
fprintf('  %-9s | peak|theta| peak|z| peak|delta| min|GM| min|PM|\n', 'controller');
for k = 1:3
    [g, p] = deal(NaN);
    if k == 2, g = min(abs(gmT)); p = min(abs(pmT)); end
    if k == 3, g = min(abs(gmQ)); p = min(abs(pmQ)); end
    fprintf('  %-9s |  %6.3f deg %6.2f m %7.3f deg  %5.2f dB %5.1f deg\n', ...
        runs(k).name, runs(k).r.pk_theta, runs(k).r.pk_z, runs(k).r.pk_delta, g, p);
end

%% ---------------------------------------------------------------- Figures
% f1 — gain-vs-q hysteresis loop
f1 = figure('Name', 'q_hysteresis', 'Color', 'w', 'Position', [80 80 760 360]);
tl1 = tiledlayout(f1, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl1, 'Gain schedule keyed on q — ascending vs descending branch');
nexttile; hold on;
plot(Qs(1:ipk), S.Kp_sched(1:ipk), 'b.-', 'LineWidth', 1.3, 'DisplayName', 'q rising (t\leq65 s)');
plot(Qs(ipk:end), S.Kp_sched(ipk:end), 'r.-', 'LineWidth', 1.3, 'DisplayName', 'q falling (t\geq65 s)');
grid on; xlabel('q [kPa]'); ylabel('K_{p,\theta}'); legend('Location', 'northwest');
title(sprintf('K_p hysteresis (\\Delta up to %.0f%%)', hyst));
nexttile; hold on;
plot(Qs(1:ipk), S.Kd_sched(1:ipk), 'b.-', 'LineWidth', 1.3, 'DisplayName', 'q rising');
plot(Qs(ipk:end), S.Kd_sched(ipk:end), 'r.-', 'LineWidth', 1.3, 'DisplayName', 'q falling');
grid on; xlabel('q [kPa]'); ylabel('K_{d,\theta}'); legend('Location', 'northwest');
title('K_d hysteresis');

% f2 — response: t-scheduled vs q-scheduled
f2 = figure('Name', 'q_response', 'Color', 'w', 'Position', [80 80 980 320]);
tl2 = tiledlayout(f2, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl2, 'Response — t-scheduled vs q-scheduled gains');
sig = {'theta', '\theta [deg]', 180/pi; 'z', 'z [m]', 1; 'delta', '\delta [deg]', 180/pi};
for j = 1:3
    nexttile; plot(runs(2).r.t, runs(2).r.(sig{j,1})*sig{j,3}, 'b-', ...
                   runs(3).r.t, runs(3).r.(sig{j,1})*sig{j,3}, 'r-', 'LineWidth', 1.2);
    xline(72, 'k:'); grid on; xlabel('t [s]'); ylabel(sig{j,2});
    if j == 1, legend('t-sched', 'q-sched', 'Location', 'best'); end
end

% f3 — frozen-time margins: t-scheduled vs q-scheduled
f3 = figure('Name', 'q_margins', 'Color', 'w', 'Position', [80 80 760 540]);
tl3 = tiledlayout(f3, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl3, 'Frozen-time margins — t-scheduled vs q-scheduled');
nexttile; plot(tm, abs(gmT), 'b.-', tm, abs(gmQ), 'r.-', 'LineWidth', 1.2);
yline(6, 'k--', 'target 6 dB'); xline(72, 'k:'); grid on;
xlabel('t [s]'); ylabel('|GM| [dB]'); legend('t-sched', 'q-sched', 'Location', 'best');
title('Gain margin');
nexttile; plot(tm, abs(pmT), 'b.-', tm, abs(pmQ), 'r.-', 'LineWidth', 1.2);
yline(30, 'k--', 'target 30^\circ'); xline(72, 'k:'); grid on;
xlabel('t [s]'); ylabel('|PM| [deg]'); legend('t-sched', 'q-sched', 'Location', 'best');
title('Phase margin');

%% Export
fig_dir = fullfile(here, 'figures');
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end
for f = [f1 f2 f3]
    try theme(f, 'light'); catch, set(f, 'Color', 'w'); end
    exportgraphics(f, fullfile(fig_dir, ['qsched_' get(f, 'Name') '.png']), 'Resolution', 200);
end
fprintf('\nFigures written to %s\n', fig_dir);

%% ------------------------------------------------------------ local helpers
function M = make_model(S, mode, fKp_q, fKd_q)
%MAKE_MODEL  ODE_LPV_ASCENT struct; mode = frozen | t (time) | q (dyn. pressure).
M = struct('fc1', S.fc1, 'fc2', S.fc2, 'fc3', S.fc3, 'fc4', S.fc4, ...
           'fc5', S.fc5, 'fc6', S.fc6, 'fc7', S.fc7, 'windfun', S.windfun, ...
           'Kp_th0', S.K0.Kp_th, 'Kd_th0', S.K0.Kd_th, ...
           'Kp_z', S.K0.Kp_z, 'Kd_z', S.K0.Kd_z, 'sched', true);
switch mode
    case 'frozen', M.sched = false; M.fKp = S.fKp;  M.fKd = S.fKd;
    case 't',      M.fKp = S.fKp;   M.fKd = S.fKd;
    case 'q',      M.fKp = fKp_q;   M.fKd = fKd_q;
end
end

function r = unpack(tt, x, S, M)
if M.sched, Kp = M.fKp(tt); Kd = M.fKd(tt);
else,       Kp = M.Kp_th0*ones(size(tt)); Kd = M.Kd_th0*ones(size(tt)); end
delta = -(Kp.*x(:, 3) + Kd.*x(:, 4) + S.K0.Kp_z*x(:, 1) + S.K0.Kd_z*x(:, 2));
r = struct('t', tt, 'theta', x(:, 3), 'z', x(:, 1), 'delta', delta);
r.pk_theta = max(abs(r.theta))*180/pi;
r.pk_z = max(abs(r.z));  r.pk_delta = max(abs(delta))*180/pi;
end

function [gm_db, pm_deg] = loop_margin(G, K)
[L, ~] = assemble_loop(G, K, []);
[g, p] = margin(L);
gm_db = 20*log10(g);  pm_deg = p;
end
