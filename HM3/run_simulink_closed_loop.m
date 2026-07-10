function out = run_simulink_closed_loop(task, o)
% Simulate hm3_closed_loop.slx and overlay it on the script response.
%   INPUT
%     task - 1 | 2 | 3 (see init_simulink_hm3)
%     o    - name-value mu_alpha_scale/mu_c_scale/severity/profile, forwarded
%            to init_simulink_hm3
%   OUTPUT
%     out - struct('script', rs, 'simulink', sl); [] if the model is missing
%   Writes figures/task<task>_simulink_vs_script.png. The .slx is built by
%   hand (models/SIMULINK_GUIDE.md), not auto-generated; if absent, prints the
%   next step and returns []. Model must log: theta_sl z_sl zdot_sl delta_sl.

arguments
    task (1,1) {mustBeMember(task, [1 2 3])} = 2
    o.mu_alpha_scale (1,1) {mustBeNumeric, mustBeReal} = 1.0
    o.mu_c_scale     (1,1) {mustBeNumeric, mustBeReal} = 1.0
    o.severity {mustBeMember(o.severity, ["light","moderate","severe"])} = 'severe'
    o.profile  {mustBeTextScalar} = 'gust'
end

here  = fileparts(mfilename('fullpath'));
model = 'hm3_closed_loop';
mdlfile = fullfile(here,'models',[model '.slx']);
if ~isfile(mdlfile)
    fprintf(['[run_simulink_closed_loop] Model not found:\n  %s\n' ...
             'Build it first by following models/SIMULINK_GUIDE.md, then re-run.\n'], mdlfile);
    out = [];
    return;
end

%% Init workspace (gains, matrices, wind, Tstop) and the script baseline
optArgs = namedargs2cell(o);
S = init_simulink_hm3(task, optArgs{:});

% script baseline (same controller/plant the model embeds)
p = S.p;
if task == 1
    G = build_plant_rigid(p);  Wact = [];
else
    G = build_plant_full(p,'ins');
    Wact = build_tvc(p,3) * build_notch_filter(p.wBM,0.002,0.7,+1);
end
K.Kp_th=S.Kp_th; K.Kd_th=S.Kd_th; K.Kp_z=S.Kp_z; K.Kd_z=S.Kd_z;
[~,T] = assemble_loop(G,K,Wact);
% replay the exact wind the model sees (whatever profile/severity was set)
w = struct('t', S.wind_ts.Time(:), 'alphaw', squeeze(S.wind_ts.Data), 'V', p.V);
rs = simulate_gust_response(T,w);

%% Simulate the Simulink model
addpath(fullfile(here,'models'));
so = sim(model, 'StopTime', num2str(S.Tstop));

get_ts = @(nm) get_logged_signal(so, nm);
[sl.t, sl.theta] = get_ts('theta_sl');
[~,    sl.z]     = get_ts('z_sl');
[~,    sl.zdot]  = get_ts('zdot_sl');
[~,    sl.delta] = get_ts('delta_sl');

%% Overlay
f = figure('Name','simulink_vs_script','Color','w','Position',[100 100 820 360]);
tl = tiledlayout(f,1,3,'TileSpacing','compact','Padding','compact');
title(tl,sprintf('Task %d - Simulink vs script (closed loop)',task));
nexttile; plot(rs.t,rs.theta*180/pi,'b-',sl.t,sl.theta*180/pi,'r--','LineWidth',1.3);
grid on; xlabel('t [s]'); ylabel('\theta [deg]'); legend('script','Simulink');
nexttile; plot(rs.t,rs.z,'b-',sl.t,sl.z,'r--','LineWidth',1.3);
grid on; xlabel('t [s]'); ylabel('z [m]'); legend('script','Simulink');
nexttile; plot(rs.t,rs.delta*180/pi,'b-',sl.t,sl.delta*180/pi,'r--','LineWidth',1.3);
grid on; xlabel('t [s]'); ylabel('\delta [deg]'); legend('script','Simulink');

fig_dir = fullfile(here,'figures');
if ~exist(fig_dir,'dir'); mkdir(fig_dir); end
try
    theme(f, 'light');    % force light theme (ignore desktop dark mode)
catch
    f.Color = 'w';        % fallback for pre-R2025a MATLAB
end
suffix = '';                       % non-default wind -> separate figure file
if ~strcmpi(o.profile, 'gust')     % keep the default gust figure unsuffixed
    suffix = ['_' lower(char(o.profile))];
end
exportgraphics(f, fullfile(fig_dir,sprintf('task%d_simulink_vs_script%s.png',task,suffix)),'Resolution',200);

out = struct('script',rs,'simulink',sl);
end

function [t,y] = get_logged_signal(so, name)
% Fetch a logged signal from a SimulationOutput by name.
%   INPUT  so, name - SimulationOutput, signal name
%   OUTPUT t, y      - time/data ([] if not found)
t = []; y = [];
try
    if isprop(so,'logsout') && ~isempty(so.logsout) && ...
            any(strcmp(so.logsout.getElementNames, name))
        e = so.logsout.getElement(name);
        t = e.Values.Time;  y = e.Values.Data;
    elseif isprop(so,name) || isfield(so,name)
        v = so.(name);  t = v.Time;  y = v.Data;
    end
catch
    warning('run_simulink_closed_loop:signal','Could not read signal "%s".',name);
end
end
