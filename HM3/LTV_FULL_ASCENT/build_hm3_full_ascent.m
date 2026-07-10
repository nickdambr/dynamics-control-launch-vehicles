function build_hm3_full_ascent(o)
% Author hm3_full_ascent.slx from code: time-varying (LPV) rigid pitch-plane
% LV, with the professor's wind generator in the loop. Simulink mirror of
% ODE_LPV_ASCENT / MAIN_FULL_ASCENT (source of truth). Run INIT_SIMULINK_LPV
% first so the referenced base variables exist.
%   INPUT
%     o.open - open model after build (default false)
%
% Layout: Clock -> 1-D lookups c1..c7, invV (bp lpv_t); plant as products +
% integrators; WindGen = one-time copy of strong_wind/Subsystem (read-only,
% never saved); controller delta = -(Kp*theta + Kd*thetadot + Kp_z*z + Kd_z*zdot),
% pitch gains switched FROZEN<->SCHEDULE by scalar 'sched' (0/1).
%
%   See also INIT_SIMULINK_LPV, RUN_FULL_ASCENT_SIMULINK, ODE_LPV_ASCENT.

arguments
    o.open (1,1) logical = false
end

mdl  = 'hm3_full_ascent';
here = fileparts(mfilename('fullpath'));
hm3  = fileparts(here);
gdir = fullfile(hm3, 'General', 'hw3-v3');
addpath(here); addpath(hm3);

%% Load the professor's generator (read-only) to copy its Subsystem
load_system(fullfile(gdir, 'strong_wind.slx'));
cleanupSW = onCleanup(@() close_system('strong_wind', 0));      % never saved

%% Fresh model + solver
if bdIsLoaded(mdl), close_system(mdl, 0); end
new_system(mdl);
% MaxStep bounded to resolve the generator's 0.1 s noise innovations;
% ode45 replay then overlays to ~1e-7 rad on theta.
set_param(mdl, 'Solver', 'ode45', 'SolverType', 'Variable-step', ...
          'StartTime', '0', 'StopTime', 'Tstop', 'MaxStep', '0.02', ...
          'SignalLogging', 'on', 'SignalLoggingName', 'logsout');

% block-add / wire helpers
A = @(name, src, pos, varargin) add_block(src, [mdl '/' name], 'Position', pos, varargin{:});
W = @(s, d) add_line(mdl, s, d, 'autorouting', 'on');
LK = 'simulink/Lookup Tables/n-D Lookup Table';

%% Clock
A('Clock', 'built-in/Clock', [30 400 60 430]);

%% Coefficient lookups on flight time (c1..c7, invV)
cs = {'c1','c2','c3','c4','c5','c6','c7','invV'};
tb = {'lpv_c1','lpv_c2','lpv_c3','lpv_c4','lpv_c5','lpv_c6','lpv_c7','lpv_invV'};
for k = 1:numel(cs)
    A(cs{k}, LK, [140 40+40*k 190 70+40*k], ...
      'NumberOfTableDimensions','1','BreakpointsForDimension1','lpv_t','Table',tb{k});
    W('Clock/1', [cs{k} '/1']);
end

%% Wind generator (one-time copy) + alpha_w = (v_wp + turb)*invV
add_block('strong_wind/Subsystem', [mdl '/WindGen'], 'Position', [140 470 220 540]);
A('vw_sum',  'built-in/Sum',     [260 480 290 540], 'Inputs','++','IconShape','rectangular');
A('aw_prod', 'built-in/Product', [330 490 360 530], 'Inputs','2');
W('Clock/1', 'WindGen/1');
W('WindGen/1', 'vw_sum/1');
W('WindGen/2', 'vw_sum/2');
W('vw_sum/1',  'aw_prod/1');
W('invV/1',    'aw_prod/2');

%% Plant: products  ->  sums  ->  integrators
% zddot products
A('P1','built-in/Product',[430 60 460 90],'Inputs','2');   % c1*zdot
A('P2','built-in/Product',[430 100 460 130],'Inputs','2'); % c2*theta
A('P3','built-in/Product',[430 140 460 170],'Inputs','2'); % c3*delta
A('P4','built-in/Product',[430 180 460 210],'Inputs','2'); % c4*alpha_w
A('zdd','built-in/Sum',[520 110 550 200],'Inputs','+++-','IconShape','rectangular');
A('int_zd','built-in/Integrator',[590 130 620 160],'InitialCondition','0'); % -> zdot
A('int_z', 'built-in/Integrator',[660 130 690 160],'InitialCondition','0'); % -> z
% thetaddot products
A('P5','built-in/Product',[430 250 460 280],'Inputs','2'); % c5*zdot
A('P6','built-in/Product',[430 290 460 320],'Inputs','2'); % c6*theta
A('P7','built-in/Product',[430 330 460 360],'Inputs','2'); % c7*delta
A('P8','built-in/Product',[430 370 460 400],'Inputs','2'); % c6*alpha_w
A('thdd','built-in/Sum',[520 300 550 390],'Inputs','+++-','IconShape','rectangular');
A('int_thd','built-in/Integrator',[590 320 620 350],'InitialCondition','0'); % -> thetadot
A('int_th', 'built-in/Integrator',[660 320 690 350],'InitialCondition','0'); % -> theta

% coefficient -> product (port 1)
W('c1/1','P1/1'); W('c2/1','P2/1'); W('c3/1','P3/1'); W('c4/1','P4/1');
W('c5/1','P5/1'); W('c6/1','P6/1'); W('c6/1','P8/1'); W('c7/1','P7/1');
% products -> sums
W('P1/1','zdd/1'); W('P2/1','zdd/2'); W('P3/1','zdd/3'); W('P4/1','zdd/4');
W('P5/1','thdd/1'); W('P6/1','thdd/2'); W('P7/1','thdd/3'); W('P8/1','thdd/4');
% integrator chains
W('zdd/1','int_zd/1');  W('int_zd/1','int_z/1');
W('thdd/1','int_thd/1'); W('int_thd/1','int_th/1');

% state feedback into the products (port 2)
W('int_zd/1','P1/2');   % zdot * c1
W('int_zd/1','P5/2');   % zdot * c5
W('int_th/1','P2/2');   % theta * c2
W('int_th/1','P6/2');   % theta * c6
W('aw_prod/1','P4/2');  % alpha_w * c4
W('aw_prod/1','P8/2');  % alpha_w * c6

%% Controller: delta = -(Kp*theta + Kd*thetadot + Kp_z*z + Kd_z*zdot)
A('sched','built-in/Constant',[760 600 800 630],'Value','sched');
A('Kp_f','built-in/Constant',[760 470 800 500],'Value','Kp_th0');
A('Kd_f','built-in/Constant',[760 540 800 570],'Value','Kd_th0');
% Clip below tsched(1): hold endpoint gain, matching the ode replay's
% 'nearest' extrapolation so the scheduled overlay is exact.
A('Kp_s',LK,[760 410 810 440],'NumberOfTableDimensions','1','BreakpointsForDimension1','tsched','Table','Kp_sched','ExtrapMethod','Clip');
A('Kd_s',LK,[760 660 810 690],'NumberOfTableDimensions','1','BreakpointsForDimension1','tsched','Table','Kd_sched','ExtrapMethod','Clip');
% Switch: u2 >= 0.5 -> port1 (scheduled), else port3 (frozen)
A('Kp_sw','built-in/Switch',[870 410 900 480],'Criteria','u2 >= Threshold','Threshold','0.5');
A('Kd_sw','built-in/Switch',[870 540 900 610],'Criteria','u2 >= Threshold','Threshold','0.5');
W('Kp_s/1','Kp_sw/1'); W('sched/1','Kp_sw/2'); W('Kp_f/1','Kp_sw/3');
W('Kd_s/1','Kd_sw/1'); W('sched/1','Kd_sw/2'); W('Kd_f/1','Kd_sw/3');
W('Clock/1','Kp_s/1'); W('Clock/1','Kd_s/1');

A('Kp_prod','built-in/Product',[940 420 970 450],'Inputs','2');  % Kp*theta
A('Kd_prod','built-in/Product',[940 550 970 580],'Inputs','2');  % Kd*thetadot
A('Gain_Kpz','built-in/Gain',[940 640 980 670],'Gain','Kp_z0');  % Kp_z*z
A('Gain_Kdz','built-in/Gain',[940 700 980 730],'Gain','Kd_z0');  % Kd_z*zdot
A('delta_sum','built-in/Sum',[1030 500 1060 620],'Inputs','----','IconShape','rectangular');
W('Kp_sw/1','Kp_prod/1'); W('int_th/1','Kp_prod/2');
W('Kd_sw/1','Kd_prod/1'); W('int_thd/1','Kd_prod/2');
W('int_z/1','Gain_Kpz/1');
W('int_zd/1','Gain_Kdz/1');
W('Kp_prod/1','delta_sum/1');
W('Kd_prod/1','delta_sum/2');
W('Gain_Kpz/1','delta_sum/3');
W('Gain_Kdz/1','delta_sum/4');

% delta back into the plant (port 2 of c3*delta, c7*delta)
W('delta_sum/1','P3/2');
W('delta_sum/1','P7/2');

%% Logging (contract with RUN_FULL_ASCENT_SIMULINK)
tw = {'-1','MaxDataPoints','inf'};   % every solver step, no 1000-pt cap
A('log_theta','built-in/ToWorkspace',[760 110 820 140],'VariableName','theta_sl','SaveFormat','Timeseries','SampleTime',tw{:});
A('log_z',    'built-in/ToWorkspace',[760 160 820 190],'VariableName','z_sl','SaveFormat','Timeseries','SampleTime',tw{:});
A('log_zdot', 'built-in/ToWorkspace',[760 210 820 240],'VariableName','zdot_sl','SaveFormat','Timeseries','SampleTime',tw{:});
A('log_delta','built-in/ToWorkspace',[1120 540 1180 570],'VariableName','delta_sl','SaveFormat','Timeseries','SampleTime',tw{:});
A('log_aw',   'built-in/ToWorkspace',[430 470 490 500],'VariableName','alpha_w_sl','SaveFormat','Timeseries','SampleTime',tw{:});
W('int_th/1','log_theta/1');
W('int_z/1', 'log_z/1');
W('int_zd/1','log_zdot/1');
W('delta_sum/1','log_delta/1');
W('aw_prod/1','log_aw/1');

%% Save
mdlfile = fullfile(here, [mdl '.slx']);
save_system(mdl, mdlfile);
fprintf('build_hm3_full_ascent: wrote %s\n', mdlfile);
if o.open, open_system(mdl); else, close_system(mdl, 0); end
end
