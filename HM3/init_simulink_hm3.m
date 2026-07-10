function S = init_simulink_hm3(task, o)
% Populate the base workspace for models/hm3_closed_loop.slx.
%   INPUT
%     task - 1 rigid+ideal actuator | 2 full+TVC+delay+notch (default)
%            | 3 full, mu_alpha/mu_c corners
%     o    - name-value:
%              mu_alpha_scale, mu_c_scale  corner scaling (Task 3), default 1
%              severity  wind severity, default 'severe'
%              profile   wind profile, default 'gust' ('strongwind' = generator)
%              push      push vars to base workspace, default true
%   OUTPUT
%     S - struct of the same data. Exported vars (task-dependent):
%           A_rigid/A_full, Bdelta_* Bwind_* (delta/alpha_w columns),
%           C_meas_* C_plot_*, Kp_th Kd_th Kp_z Kd_z,
%           tvc_num/tvc_den, notch_num/notch_den, wind_ts, Tstop
%   Gains/filters are exactly those designed by the scripts; the model only
%   mirrors them. Re-run whenever a gain changes.

arguments
    task (1,1) {mustBeMember(task, [1 2 3])} = 2
    o.mu_alpha_scale (1,1) {mustBeNumeric, mustBeReal} = 1.0
    o.mu_c_scale     (1,1) {mustBeNumeric, mustBeReal} = 1.0
    o.severity {mustBeMember(o.severity, ["light","moderate","severe"])} = 'severe'
    o.profile  {mustBeTextScalar} = 'gust'
    o.push (1,1) logical = true
end

%% Parameters and controller (designed by the scripts, reused here)
%  Controller and notch FROZEN at the nominal point, as in main_task3:
%  robustness = fixed nominal gains on the perturbed plant. Re-tuning on the
%  corner plant would destabilise the bending-augmented loop.
p  = load_hw3_params('mu_alpha_scale',o.mu_alpha_scale,'mu_c_scale',o.mu_c_scale);
p0 = load_hw3_params();                                   % nominal design point
K  = design_controller(build_plant_rigid(p0), [], 'verbose', false);

S = struct();
S.p = p;
S.Kp_th = K.Kp_th; S.Kd_th = K.Kd_th; S.Kp_z = K.Kp_z; S.Kd_z = K.Kd_z;

%% Plant matrices (split input columns so Simulink can wire delta / alpha_w)
Gr = build_plant_rigid(p);
S.A_rigid  = Gr.A;
S.Bdelta_rigid = Gr.B(:,1);  S.Bwind_rigid = Gr.B(:,2);
S.C_meas_rigid = Gr.C(1:4,:);            % [theta_m thetadot_m z_m zdot_m]
S.C_plot_rigid = Gr.C(5:7,:);            % [theta z zdot]

Gf = build_plant_full(p,'ins');
S.A_full  = Gf.A;
S.Bdelta_full = Gf.B(:,1);  S.Bwind_full = Gf.B(:,2);
S.C_meas_full = Gf.C(1:4,:);
S.C_plot_full = Gf.C(5:7,:);

%% Actuator (TVC + delay) and bending notch
Wtvc = build_tvc(p,3);
[S.tvc_num, S.tvc_den] = tfdata(tf(Wtvc),'v');
Hx = build_notch_filter(p0.wBM, 0.002, 0.7, +1);
[S.notch_num, S.notch_den] = tfdata(Hx,'v');

%% Wind disturbance as a timeseries for a From Workspace block
w = load_wind_profile(p,'severity',o.severity,'profile',o.profile);
S.wind_ts = timeseries(w.alphaw(:), w.t(:), 'Name','alpha_w');
S.Tstop = w.t(end);

S.task = task;

%% Push to base workspace
if o.push
    fn = fieldnames(S);
    for i = 1:numel(fn)
        assignin('base', fn{i}, S.(fn{i}));
    end
    fprintf('init_simulink_hm3: pushed %d variables to base (task %d, mu_a=%.2f mu_c=%.2f).\n', ...
            numel(fn), task, o.mu_alpha_scale, o.mu_c_scale);
end
end
