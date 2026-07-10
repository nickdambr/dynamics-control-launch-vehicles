function out = run_flex_simulink(o)
% Simulate hm3_full_ascent_flex.slx (bending + INS + TVC + varying notch,
% frozen gains) and overlay vs the ode45 baseline (ODE_LPV_FLEX), driven by
% the wind the model generated. Writes flex_simulink_vs_script.png.
%   INPUT
%     o.rebuild - re-author model with BUILD_HM3_FULL_ASCENT_FLEX first (default false)
%   OUTPUT
%     out - struct: t, err (theta/eta/delta)
%
%   See also BUILD_HM3_FULL_ASCENT_FLEX, ODE_LPV_FLEX, MAIN_FLEX.

arguments
    o.rebuild (1,1) logical = false
end

here = fileparts(mfilename('fullpath'));
addpath(here); addpath(fileparts(here));
mdl     = 'hm3_full_ascent_flex';
mdlfile = fullfile(here, [mdl '.slx']);
if o.rebuild || ~isfile(mdlfile)
    build_hm3_full_ascent_flex();
end

S    = init_simulink_lpv();
nt   = size(S.tvc.At, 1);
odeo = odeset('RelTol', 1e-9, 'AbsTol', 1e-11);
load_system(mdlfile);

so = sim(mdl, 'StopTime', num2str(S.Tstop));
th = so.theta_sl;  et = so.eta_sl;  de = so.delta_sl;  aw = so.alpha_w_sl;
tt = th.Time;

% ode45 replay (varying notch, frozen gains) on the SAME wind
M = struct('fa1', S.fa1, 'fa3', S.fa3, 'fa4', S.fa4, 'fA6', S.fA6, 'fK1', S.fK1, ...
           'fV', S.fV, 'fomega', S.fomega, 'faqk', S.faqk, 'fsig', S.fsig, 'fphi', S.fphi, ...
           'windfun', griddedInterpolant(aw.Time, squeeze(aw.Data), 'linear', 'nearest'), ...
           'fwn', S.fomega, 'zN', S.notch.zN, 'zD', S.notch.zD, 'zBM', S.notch.zBM, ...
           'At', S.tvc.At, 'Bt', S.tvc.Bt, 'Ct', S.tvc.Ct, 'Dt', S.tvc.Dt, ...
           'fKp', S.fKp, 'fKd', S.fKd, 'Kp_th0', S.K0.Kp_th, 'Kd_th0', S.K0.Kd_th, ...
           'Kp_z', S.K0.Kp_z, 'Kd_z', S.K0.Kd_z, 'sched', false);
[~, x] = ode45(@(t, x) ode_lpv_flex(t, x, M), tt, zeros(6 + 2 + nt, 1), odeo);
delta_ode = M.Ct*x(:, 9:end).' + M.Dt*0;            % delta = Ct*x_tvc (Dt=0)
delta_ode = delta_ode(:);

err = struct('theta', max(abs(x(:, 3) - squeeze(th.Data))), ...
             'eta',   max(abs(x(:, 5) - squeeze(et.Data))), ...
             'delta', max(abs(delta_ode - squeeze(de.Data))));
out = struct('t', tt, 'err', err);
fprintf('flex overlay: max|d theta|=%.2e rad  |d eta|=%.2e  |d delta|=%.2e rad\n', ...
        err.theta, err.eta, err.delta);

%% Overlay figure (baseline solid, Simulink dashed)
f = figure('Name', 'flex_simulink_vs_script', 'Color', 'w', 'Position', [80 80 1000 320]);
tl = tiledlayout(f, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Flexible LPV — Simulink vs ode45 (varying notch)   max|\\Delta\\theta| = %.1e rad', err.theta));
nexttile; plot(tt, x(:, 3)*180/pi, 'b-', th.Time, squeeze(th.Data)*180/pi, 'r--', 'LineWidth', 1.1);
grid on; xlabel('t [s]'); ylabel('\theta [deg]'); legend('ode45', 'Simulink', 'Location', 'best');
nexttile; plot(tt, x(:, 5), 'b-', et.Time, squeeze(et.Data), 'r--', 'LineWidth', 1.1);
grid on; xlabel('t [s]'); ylabel('\eta (bending)'); title('bending coordinate');
nexttile; plot(tt, delta_ode*180/pi, 'b-', de.Time, squeeze(de.Data)*180/pi, 'r--', 'LineWidth', 1.1);
grid on; xlabel('t [s]'); ylabel('\delta [deg]');

fig_dir = fullfile(here, 'figures');
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end
try theme(f, 'light'); catch, set(f, 'Color', 'w'); end
exportgraphics(f, fullfile(fig_dir, 'flex_simulink_vs_script.png'), 'Resolution', 200);
fprintf('Overlay figure written to %s\n', fig_dir);
end
