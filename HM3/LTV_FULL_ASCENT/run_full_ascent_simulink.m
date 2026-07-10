function out = run_full_ascent_simulink(o)
% Simulate hm3_full_ascent.slx (frozen + scheduled) and overlay vs the ode45
% baseline (ODE_LPV_ASCENT). Replay driven by the wind the model generated
% (logged alpha_w), so the residual is just wind interpolation between solver
% steps (~1e-7 rad on theta). Writes fullascent_simulink_vs_script.png.
% Model logs (To Workspace, Timeseries): theta_sl, z_sl, zdot_sl, delta_sl,
% alpha_w_sl.
%   INPUT
%     o.rebuild - re-author model with BUILD_HM3_FULL_ASCENT first (default false)
%   OUTPUT
%     out - struct: frozen/scheduled, each with t and err (theta/z/delta)
%
%   See also INIT_SIMULINK_LPV, BUILD_HM3_FULL_ASCENT, ODE_LPV_ASCENT,
%   MAIN_FULL_ASCENT.

arguments
    o.rebuild (1,1) logical = false
end

here = fileparts(mfilename('fullpath'));
addpath(here); addpath(fileparts(here));
mdl     = 'hm3_full_ascent';
mdlfile = fullfile(here, [mdl '.slx']);
if o.rebuild || ~isfile(mdlfile)
    build_hm3_full_ascent();
end

S    = init_simulink_lpv();
odeo = odeset('RelTol', 1e-9, 'AbsTol', 1e-11);
load_system(mdlfile);

variants = {'frozen', 0; 'scheduled', 1};
out = struct();

f = figure('Name', 'simulink_vs_script', 'Color', 'w', 'Position', [80 80 1000 560]);
tl = tiledlayout(f, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'Full-ascent LPV — Simulink vs MATLAB ode45 baseline');

for v = 1:2
    nm = variants{v, 1};  sc = variants{v, 2};
    assignin('base', 'sched', sc);
    so = sim(mdl, 'StopTime', num2str(S.Tstop));

    th = so.theta_sl;  z = so.z_sl;  de = so.delta_sl;  aw = so.alpha_w_sl;
    tt = th.Time;

    % ode45 replay on the SAME wind the model produced
    M = pack_model(S, sc, aw);
    [~, x] = ode45(@(t, x) ode_lpv_ascent(t, x, M), tt, zeros(4, 1), odeo);
    if sc, Kp = S.fKp(tt); Kd = S.fKd(tt);
    else,  Kp = S.K0.Kp_th*ones(size(tt)); Kd = S.K0.Kd_th*ones(size(tt)); end
    del_ode = -(Kp.*x(:, 3) + Kd.*x(:, 4) + S.K0.Kp_z*x(:, 1) + S.K0.Kd_z*x(:, 2));

    err = struct('theta', max(abs(x(:, 3) - squeeze(th.Data))), ...
                 'z',     max(abs(x(:, 1) - squeeze(z.Data))), ...
                 'delta', max(abs(del_ode - squeeze(de.Data))));
    out.(nm) = struct('t', tt, 'err', err);
    fprintf('%-9s overlay: max|d theta|=%.2e rad  |d z|=%.2e m  |d delta|=%.2e rad\n', ...
            nm, err.theta, err.z, err.delta);

    % overlay (baseline solid, Simulink dashed)
    nexttile; plot(tt, x(:, 3)*180/pi, 'b-', th.Time, squeeze(th.Data)*180/pi, 'r--', 'LineWidth', 1.1);
    grid on; xlabel('t [s]'); ylabel('\theta [deg]');
    title(sprintf('%s: \\theta (\\Delta=%.1e rad)', nm, err.theta));
    if v == 1, legend('ode45', 'Simulink', 'Location', 'best'); end
    nexttile; plot(tt, x(:, 1), 'b-', th.Time, squeeze(z.Data), 'r--', 'LineWidth', 1.1);
    grid on; xlabel('t [s]'); ylabel('z [m]'); title([nm ': z']);
    nexttile; plot(tt, del_ode*180/pi, 'b-', th.Time, squeeze(de.Data)*180/pi, 'r--', 'LineWidth', 1.1);
    grid on; xlabel('t [s]'); ylabel('\delta [deg]'); title([nm ': \delta']);
end

fig_dir = fullfile(here, 'figures');
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end
try theme(f, 'light'); catch, set(f, 'Color', 'w'); end
exportgraphics(f, fullfile(fig_dir, 'fullascent_simulink_vs_script.png'), 'Resolution', 200);
fprintf('Overlay figure written to %s\n', fig_dir);
end

% ------------------------------------------------------------------------
function M = pack_model(S, sched, aw)
% Build the ODE_LPV_ASCENT struct, wind = the model's logged alpha_w.
%   INPUT
%     S     - init_simulink_lpv struct
%     sched - 0 frozen gains, 1 scheduled
%     aw    - logged alpha_w timeseries
%   OUTPUT
%     M - struct for ODE_LPV_ASCENT
M = struct('fc1', S.fc1, 'fc2', S.fc2, 'fc3', S.fc3, 'fc4', S.fc4, ...
           'fc5', S.fc5, 'fc6', S.fc6, 'fc7', S.fc7, 'fKp', S.fKp, 'fKd', S.fKd, ...
           'Kp_th0', S.K0.Kp_th, 'Kd_th0', S.K0.Kd_th, 'Kp_z', S.K0.Kp_z, ...
           'Kd_z', S.K0.Kd_z, 'sched', logical(sched), ...
           'windfun', griddedInterpolant(aw.Time, squeeze(aw.Data), 'linear', 'nearest'));
end
