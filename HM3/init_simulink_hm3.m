function S = init_simulink_hm3(task, varargin)
%INIT_SIMULINK_HM3  Populate the base workspace for the Simulink model.
%
%   S = INIT_SIMULINK_HM3(task) computes everything the block diagram
%   models/hm3_closed_loop.slx needs and pushes it to the base workspace as
%   named variables, so the State-Space / Transfer-Fcn / Gain blocks can
%   reference them. It returns the same data as a struct S for convenience.
%
%   The script is the single source of truth: the gains and filters are the
%   very ones designed and validated by the MATLAB scripts (Task 1-3), so
%   the Simulink model only mirrors them. Run this before opening/simulating
%   the model, and re-run it whenever a gain changes.
%
%   task = 1  rigid plant + ideal actuator (Task 1)
%        = 2  full plant + TVC + delay + notch (Task 2)         [default]
%        = 3  full plant, parametrised by mu_alpha/mu_c corners (Task 3)
%
%   Name/value options:
%     'mu_alpha_scale', 'mu_c_scale'   corner scaling (Task 3), default 1
%     'severity'                        wind severity, default 'severe'
%     'profile'                         wind profile, default 'gust';
%                                       'strongwind' = professor's generator
%     'push', true/false                push to base workspace, default true
%
%   Variables exported (subset, depending on task):
%     A_rigid B_rigid C_meas C_plot          rigid plant (Task 1)
%     A_full  B_full                          full 6-state plant (Task 2/3)
%     Bdelta_* Bwind_*                        input columns (delta / alpha_w)
%     Kp_th Kd_th Kp_z Kd_z                   controller gains
%     tvc_num tvc_den                         TVC + delay transfer function
%     notch_num notch_den                     bending notch (Eq. 4)
%     wind_ts                                 alpha_w(t) timeseries (From Wkspc)
%     Tstop                                   suggested simulation stop time
%
%   See also LOAD_HW3_PARAMS, BUILD_PLANT_RIGID, BUILD_PLANT_FULL,
%   DESIGN_CONTROLLER, SIMULINK_GUIDE (models/SIMULINK_GUIDE.md).

if nargin < 1 || isempty(task), task = 2; end
ip = inputParser;
ip.addParameter('mu_alpha_scale',1.0);
ip.addParameter('mu_c_scale',1.0);
ip.addParameter('severity','severe');
ip.addParameter('profile','gust');
ip.addParameter('push',true);
ip.parse(varargin{:});
o = ip.Results;

%% Parameters and controller (designed by the scripts, reused here)
%  The controller and the notch are FROZEN at the nominal design point, as
%  in main_task3: robustness means fixed nominal gains evaluated on the
%  perturbed plant. Re-tuning on the corner plant would destabilize the
%  bending-augmented loop (and is not what the assignment asks).
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
