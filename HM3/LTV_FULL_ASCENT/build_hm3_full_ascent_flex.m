function build_hm3_full_ascent_flex(o)
% Author hm3_full_ascent_flex.slx: BUILD_HM3_FULL_ASCENT plus the flexible
% architecture (first bending mode, INS coupling, TVC + transport delay,
% VARYING notch tracking omega(t)). Simulink mirror of ODE_LPV_FLEX.
% PD gains FROZEN at max-q; the showcase is notch tracking, not scheduling.
% Run INIT_SIMULINK_LPV first (pushes lpv_omega, lpv_2zBMw, lpv_aqk, lpv_sig,
% lpv_phi, notch_zN/zD, tvc_num/tvc_den, ...).
%   INPUT
%     o.open - open model after build (default false)
%
% z/theta plant identical to the rigid model (bending does not feed it).
% Added: bending etaddot = -omega^2*eta - 2*zBM*omega*etadot + aqk*delta;
% INS theta_m = theta + sig*eta, z_m = z - phi*eta (and rates);
% varying notch (controllable-canonical, centre omega(t)); TVC = LTI Transfer Fcn.
%
%   See also BUILD_HM3_FULL_ASCENT, ODE_LPV_FLEX, RUN_FLEX_SIMULINK.

arguments
    o.open (1,1) logical = false
end

mdl  = 'hm3_full_ascent_flex';
here = fileparts(mfilename('fullpath'));
hm3  = fileparts(here);
gdir = fullfile(hm3, 'General', 'hw3-v3');
addpath(here); addpath(hm3);

load_system(fullfile(gdir, 'strong_wind.slx'));               % read-only source
cleanupSW = onCleanup(@() close_system('strong_wind', 0));

if bdIsLoaded(mdl), close_system(mdl, 0); end
new_system(mdl);
set_param(mdl, 'Solver', 'ode45', 'SolverType', 'Variable-step', ...
          'StartTime', '0', 'StopTime', 'Tstop', 'MaxStep', '0.02', ...
          'SignalLogging', 'on', 'SignalLoggingName', 'logsout');

A  = @(name, src, pos, varargin) add_block(src, [mdl '/' name], 'Position', pos, varargin{:});
W  = @(s, d) add_line(mdl, s, d, 'autorouting', 'on');
LK = 'simulink/Lookup Tables/n-D Lookup Table';
lk = @(name, tbl, bp, pos) A(name, LK, pos, 'NumberOfTableDimensions','1', ...
                            'BreakpointsForDimension1', bp, 'Table', tbl, 'ExtrapMethod','Clip');

%% Clock + coefficient lookups (z/theta plant, bending, INS, notch)
A('Clock', 'built-in/Clock', [20 20 50 50]);
names = {'c1','c2','c3','c4','c5','c6','c7','invV','omega2','w2zBM','aqk','sig','phi','omega'};
tbls  = {'lpv_c1','lpv_c2','lpv_c3','lpv_c4','lpv_c5','lpv_c6','lpv_c7','lpv_invV', ...
         'lpv_omega2','lpv_2zBMw','lpv_aqk','lpv_sig','lpv_phi','lpv_omega'};
for k = 1:numel(names)
    lk(names{k}, tbls{k}, 'lpv_t', [120 20+34*k 165 44+34*k]);
    W('Clock/1', [names{k} '/1']);
end

%% Wind generator (copy) -> alpha_w = (v_wp + turb)*invV
add_block('strong_wind/Subsystem', [mdl '/WindGen'], 'Position', [120 560 200 630]);
A('vw_sum',  'built-in/Sum',     [240 570 270 630], 'Inputs','++','IconShape','rectangular');
A('aw_prod', 'built-in/Product', [310 580 340 620], 'Inputs','2');
W('Clock/1','WindGen/1'); W('WindGen/1','vw_sum/1'); W('WindGen/2','vw_sum/2');
W('vw_sum/1','aw_prod/1'); W('invV/1','aw_prod/2');

%% z/theta plant (identical to the rigid model)
A('P1','built-in/Product',[430 40 460 70],'Inputs','2');
A('P2','built-in/Product',[430 80 460 110],'Inputs','2');
A('P3','built-in/Product',[430 120 460 150],'Inputs','2');
A('P4','built-in/Product',[430 160 460 190],'Inputs','2');
A('zdd','built-in/Sum',[520 90 550 180],'Inputs','+++-','IconShape','rectangular');
A('int_zd','built-in/Integrator',[590 110 620 140],'InitialCondition','0');
A('int_z', 'built-in/Integrator',[660 110 690 140],'InitialCondition','0');
A('P5','built-in/Product',[430 230 460 260],'Inputs','2');
A('P6','built-in/Product',[430 270 460 300],'Inputs','2');
A('P7','built-in/Product',[430 310 460 340],'Inputs','2');
A('P8','built-in/Product',[430 350 460 380],'Inputs','2');
A('thdd','built-in/Sum',[520 280 550 370],'Inputs','+++-','IconShape','rectangular');
A('int_thd','built-in/Integrator',[590 300 620 330],'InitialCondition','0');
A('int_th', 'built-in/Integrator',[660 300 690 330],'InitialCondition','0');
W('c1/1','P1/1'); W('c2/1','P2/1'); W('c3/1','P3/1'); W('c4/1','P4/1');
W('c5/1','P5/1'); W('c6/1','P6/1'); W('c6/1','P8/1'); W('c7/1','P7/1');
W('P1/1','zdd/1'); W('P2/1','zdd/2'); W('P3/1','zdd/3'); W('P4/1','zdd/4');
W('P5/1','thdd/1'); W('P6/1','thdd/2'); W('P7/1','thdd/3'); W('P8/1','thdd/4');
W('zdd/1','int_zd/1'); W('int_zd/1','int_z/1');
W('thdd/1','int_thd/1'); W('int_thd/1','int_th/1');
W('int_zd/1','P1/2'); W('int_zd/1','P5/2');
W('int_th/1','P2/2'); W('int_th/1','P6/2');
W('aw_prod/1','P4/2'); W('aw_prod/1','P8/2');

%% Bending: etaddot = -omega^2*eta - 2*zBM*omega*etadot + aqk*delta
A('Pe1','built-in/Product',[430 430 460 460],'Inputs','2');   % omega^2 * eta
A('Pe2','built-in/Product',[430 470 460 500],'Inputs','2');   % 2*zBM*omega * etadot
A('Pe3','built-in/Product',[430 510 460 540],'Inputs','2');   % aqk * delta
A('etadd','built-in/Sum',[520 450 550 530],'Inputs','--+','IconShape','rectangular');
A('int_etad','built-in/Integrator',[590 470 620 500],'InitialCondition','0'); % -> etadot
A('int_eta', 'built-in/Integrator',[660 470 690 500],'InitialCondition','0'); % -> eta
W('omega2/1','Pe1/1'); W('int_eta/1','Pe1/2');
W('w2zBM/1','Pe2/1');  W('int_etad/1','Pe2/2');
W('aqk/1','Pe3/1');
W('Pe1/1','etadd/1'); W('Pe2/1','etadd/2'); W('Pe3/1','etadd/3');
W('etadd/1','int_etad/1'); W('int_etad/1','int_eta/1');

%% INS measurements: theta_m=theta+sig*eta, z_m=z-phi*eta (and rates)
A('Ps1','built-in/Product',[760 410 790 440],'Inputs','2');   % sig*eta
A('Ps2','built-in/Product',[760 450 790 480],'Inputs','2');   % sig*etadot
A('Pp1','built-in/Product',[760 490 790 520],'Inputs','2');   % phi*eta
A('Pp2','built-in/Product',[760 530 790 560],'Inputs','2');   % phi*etadot
W('sig/1','Ps1/1'); W('int_eta/1','Ps1/2');
W('sig/1','Ps2/1'); W('int_etad/1','Ps2/2');
W('phi/1','Pp1/1'); W('int_eta/1','Pp1/2');
W('phi/1','Pp2/1'); W('int_etad/1','Pp2/2');
A('theta_m','built-in/Sum',[850 300 880 330],'Inputs','++','IconShape','rectangular');
A('thdot_m','built-in/Sum',[850 360 880 390],'Inputs','++','IconShape','rectangular');
A('z_m',    'built-in/Sum',[850 110 880 140],'Inputs','+-','IconShape','rectangular');
A('zdot_m', 'built-in/Sum',[850 170 880 200],'Inputs','+-','IconShape','rectangular');
W('int_th/1','theta_m/1');  W('Ps1/1','theta_m/2');
W('int_thd/1','thdot_m/1'); W('Ps2/1','thdot_m/2');
W('int_z/1','z_m/1');       W('Pp1/1','z_m/2');
W('int_zd/1','zdot_m/1');   W('Pp2/1','zdot_m/2');

%% Controller (frozen gains): u_pd = -(Kp*th_m + Kd*thd_m + Kpz*z_m + Kdz*zd_m)
A('Gth', 'built-in/Gain',[940 300 980 330],'Gain','Kp_th0');
A('Gthd','built-in/Gain',[940 360 980 390],'Gain','Kd_th0');
A('Gz',  'built-in/Gain',[940 110 980 140],'Gain','Kp_z0');
A('Gzd', 'built-in/Gain',[940 170 980 200],'Gain','Kd_z0');
A('u_pd','built-in/Sum',[1040 200 1070 320],'Inputs','----','IconShape','rectangular');
W('theta_m/1','Gth/1'); W('thdot_m/1','Gthd/1'); W('z_m/1','Gz/1'); W('zdot_m/1','Gzd/1');
W('Gz/1','u_pd/1'); W('Gzd/1','u_pd/2'); W('Gth/1','u_pd/3'); W('Gthd/1','u_pd/4');

%% Varying notch: xn1'=xn2 ; xn2'=-omega^2*xn1 - 2*zD*omega*xn2 + u_pd
A('G2zDw','built-in/Gain',[1120 40 1170 70],'Gain','2*notch_zD');          % 2*zD*omega
A('Gcout','built-in/Gain',[1120 640 1175 670],'Gain','2*(notch_zN-notch_zD)'); % output coeff
W('omega/1','G2zDw/1'); W('omega/1','Gcout/1');
A('Pn1','built-in/Product',[1200 360 1230 390],'Inputs','2');   % omega^2 * xn1
A('Pn2','built-in/Product',[1200 400 1230 430],'Inputs','2');   % 2*zD*omega * xn2
A('xn2d','built-in/Sum',[1280 360 1310 430],'Inputs','--+','IconShape','rectangular');
A('int_xn2','built-in/Integrator',[1350 380 1380 410],'InitialCondition','0'); % -> xn2
A('int_xn1','built-in/Integrator',[1420 380 1450 410],'InitialCondition','0'); % -> xn1
A('Pn3','built-in/Product',[1200 600 1230 630],'Inputs','2');   % cout * xn2
A('v_sum','built-in/Sum',[1500 300 1530 360],'Inputs','++','IconShape','rectangular');
W('omega2/1','Pn1/1'); W('int_xn1/1','Pn1/2');
W('G2zDw/1','Pn2/1');  W('int_xn2/1','Pn2/2');
W('Pn1/1','xn2d/1'); W('Pn2/1','xn2d/2'); W('u_pd/1','xn2d/3');
W('xn2d/1','int_xn2/1'); W('int_xn2/1','int_xn1/1');
W('Gcout/1','Pn3/1'); W('int_xn2/1','Pn3/2');
W('u_pd/1','v_sum/1'); W('Pn3/1','v_sum/2');

%% TVC (LTI) : v -> delta
A('TVC','built-in/TransferFcn',[1580 300 1660 360],'Numerator','tvc_num','Denominator','tvc_den');
W('v_sum/1','TVC/1');
% delta feeds plant (c3*delta, c7*delta), bending (aqk*delta), logging
W('TVC/1','P3/2'); W('TVC/1','P7/2'); W('TVC/1','Pe3/2');

%% Logging
tw = {'-1','MaxDataPoints','inf'};
A('log_theta','built-in/ToWorkspace',[760 110 820 140],'VariableName','theta_sl','SaveFormat','Timeseries','SampleTime',tw{:});
A('log_z',    'built-in/ToWorkspace',[760 160 820 190],'VariableName','z_sl','SaveFormat','Timeseries','SampleTime',tw{:});
A('log_zdot', 'built-in/ToWorkspace',[760 210 820 240],'VariableName','zdot_sl','SaveFormat','Timeseries','SampleTime',tw{:});
A('log_eta',  'built-in/ToWorkspace',[760 250 820 280],'VariableName','eta_sl','SaveFormat','Timeseries','SampleTime',tw{:});
A('log_delta','built-in/ToWorkspace',[1700 300 1760 330],'VariableName','delta_sl','SaveFormat','Timeseries','SampleTime',tw{:});
A('log_aw',   'built-in/ToWorkspace',[410 580 470 610],'VariableName','alpha_w_sl','SaveFormat','Timeseries','SampleTime',tw{:});
W('int_th/1','log_theta/1'); W('int_z/1','log_z/1'); W('int_zd/1','log_zdot/1');
W('int_eta/1','log_eta/1');  W('TVC/1','log_delta/1'); W('aw_prod/1','log_aw/1');

%% Save
mdlfile = fullfile(here, [mdl '.slx']);
save_system(mdl, mdlfile);
fprintf('build_hm3_full_ascent_flex: wrote %s\n', mdlfile);
if o.open, open_system(mdl); else, close_system(mdl, 0); end
end
