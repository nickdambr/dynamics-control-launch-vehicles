function S = init_simulink_lpv(o)
%INIT_SIMULINK_LPV  Base-workspace setup for the full-ascent LPV model.
%
%   S = INIT_SIMULINK_LPV() builds everything the full-ascent showcase needs
%   and pushes it to the base workspace as named variables, so the Simulink
%   model hm3_full_ascent.slx (lookup tables + integrators + the professor's
%   wind generator) can reference them. It also returns the same data as a
%   struct S of grid vectors and griddedInterpolant handles, ready for the
%   pure-MATLAB LTV baseline (ODE_LPV_ASCENT / MAIN_FULL_ASCENT).
%
%   This is the LPV counterpart of HM3's frozen-time INIT_SIMULINK_HM3: the
%   plant matrices of BUILD_PLANT_RIGID are no longer evaluated once at
%   t = 72 s but turned into time histories, so the wind generator and the
%   vehicle dynamics share the same clock. The rigid 4-state pitch plane is
%   used (no bending / TVC): the showcase is the gain scheduling, not the
%   flex model (see the ticket / README for the rationale).
%
%   The time-varying coefficients come from the reference data set
%   General/hw3-v3/GreensiteLPV_DATA.mat (the very file LOAD_HW3_PARAMS
%   samples at a single instant). The LPV plant, written in the same form as
%   BUILD_PLANT_RIGID,
%
%       zddot     = a1*zdot + (a1*V+a4)*theta + a3*delta - a1*V*alpha_w
%       thetaddot = (A6/V)*zdot + A6*theta    + K1*delta - A6 *alpha_w
%
%   is realised with one "effective coefficient" per term so each product is
%   a single lookup x signal (every coefficient stays inspectable):
%
%       c1 = a1        (*zdot)      c5 = A6/V   (*zdot)
%       c2 = a1*V+a4   (*theta)     c6 = A6     (*theta, *(-alpha_w))
%       c3 = a3        (*delta)     c7 = K1     (*delta)
%       c4 = a1*V      (*alpha_w)   invV = 1/V  (alpha_w = v_w*invV)
%
%   Name/value options:
%     'tsched_step'  gain-schedule grid step [s]      (default 5)
%     't0'           schedule start / sim start [s]   (default 5)
%     'Tstop'        simulation horizon [s]           (default 140)
%     'push'         push variables to base ws        (default true)
%
%   Variables exported to base (consumed by hm3_full_ascent.slx):
%     lpv_t                          coefficient breakpoints (flight time)
%     lpv_c1..lpv_c7, lpv_invV       LPV plant coefficient tables
%     lpv_Q, lpv_V, lpv_h            dyn. pressure / speed / altitude tables
%     Kp_th0 Kd_th0 Kp_z0 Kd_z0      FROZEN max-qbar PD gains (HM3 Task 1)
%     tsched Kp_sched Kd_sched       scheduled pitch gains vs flight time
%     sched                          0 = frozen, 1 = gain-scheduled (default 0)
%     drywind GreensiteLPV           inputs the wind generator needs
%     Tstart Tstop                   simulation start / stop time
%
%   See also ODE_LPV_ASCENT, MAIN_FULL_ASCENT, BUILD_HM3_FULL_ASCENT,
%   LOAD_HW3_PARAMS, BUILD_PLANT_RIGID, DESIGN_CONTROLLER.

arguments
    o.tsched_step (1,1) {mustBeNumeric, mustBeReal, mustBePositive} = 5
    o.t0          (1,1) {mustBeNumeric, mustBeReal, mustBeNonnegative} = 5
    o.Tstop       (1,1) {mustBeNumeric, mustBeReal, mustBePositive} = 140
    o.push        (1,1) logical = true
end

here = fileparts(mfilename('fullpath'));
hm3  = fileparts(here);                 % parent HM3 folder (reused helpers)
addpath(hm3);
gdir = fullfile(hm3, 'General', 'hw3-v3');

%% Reference data: LPV coefficients + dry-wind dispersion
D = load(fullfile(gdir, 'GreensiteLPV_DATA.mat'));  L = D.GreensiteLPV;
W = load(fullfile(gdir, 'drywind.mat'));            drywind = W.drywind;
GreensiteLPV = L;                                   % the generator needs it too

%% Coefficient breakpoints and effective-coefficient tables (horizon [0,Tstop])
tg   = L.V.Time(:);
tg   = tg(tg <= o.Tstop + eps);
at   = @(f) interp1(L.(f).Time, squeeze(L.(f).Data), tg);
V  = at('V');  A6 = at('A6'); K1 = at('K1');
a1 = at('a1'); a3 = at('a3'); a4 = at('a4');
Q  = at('Q');  h  = at('h');
Vsafe = max(V, 1);                      % V(0)=0: guard A6/V and 1/V at lift-off

c1   = a1;                 % * zdot
c2   = a1.*V + a4;         % * theta
c3   = a3;                 % * delta
c4   = a1.*Vsafe;          % * alpha_w  (Vsafe cancels alpha_w = v_w/Vsafe)
c5   = A6./Vsafe;          % * zdot
c6   = A6;                 % * theta and * (-alpha_w)
c7   = K1;                 % * delta
invV = 1./Vsafe;           % alpha_w = v_w * invV

%% Frozen controller: the HM3 Task-1 max-qbar design, reused unchanged
p0 = load_hw3_params();                                  % t_ref = 72 s
K0 = design_controller(build_plant_rigid(p0), [], 'verbose', false);

%% Gain schedule: design one frozen PD per grid point, warm-started (continuation)
tsched   = (o.t0:o.tsched_step:o.Tstop).';
Kp_sched = zeros(size(tsched));  Kd_sched = Kp_sched;
warnState = warning('off', 'Control:analysis:MarginUnstable');
restoreWarn = onCleanup(@() warning(warnState));
Kprev = [2.0 1.4];                                       % first warm start
for i = 1:numel(tsched)
    pk = load_hw3_params('t_ref', tsched(i));
    Kk = design_controller(build_plant_rigid(pk), [], 'K0', Kprev, 'verbose', false);
    Kp_sched(i) = Kk.Kp_th;  Kd_sched(i) = Kk.Kd_th;
    Kprev = [Kk.Kp_th Kk.Kd_th];                         % continuation
end

%% Wind: run the professor's generator ONCE over the whole ascent
wg     = run_wind_generator(gdir, o.Tstop, drywind, L);
Vwg    = max(interp1(L.V.Time, squeeze(L.V.Data), wg.t), 1);
alphaw = wg.vw ./ Vwg;                                   % wind angle of attack

%% Assemble the return struct (grids + interpolants for the LTV baseline)
S = struct();
S.tg = tg; S.V = V; S.A6 = A6; S.K1 = K1; S.a1 = a1; S.a3 = a3; S.a4 = a4;
S.Q = Q; S.h = h;
S.t0 = o.t0; S.Tstop = o.Tstop;
S.K0 = K0;                                               % frozen gains struct
S.tsched = tsched; S.Kp_sched = Kp_sched; S.Kd_sched = Kd_sched;
S.wind = struct('t', wg.t, 'vw', wg.vw, 'alphaw', alphaw);

gi = @(y) griddedInterpolant(tg, y, 'linear', 'nearest');
S.fc1 = gi(c1); S.fc2 = gi(c2); S.fc3 = gi(c3); S.fc4 = gi(c4);
S.fc5 = gi(c5); S.fc6 = gi(c6); S.fc7 = gi(c7);
S.fV  = gi(V);  S.fQ = gi(Q);
S.fKp = griddedInterpolant(tsched, Kp_sched, 'linear', 'nearest');
S.fKd = griddedInterpolant(tsched, Kd_sched, 'linear', 'nearest');
S.windfun = griddedInterpolant(wg.t, alphaw, 'linear', 'nearest');

%% Push to base workspace for the Simulink model
if o.push
    base = struct( ...
        'lpv_t', tg, 'lpv_c1', c1, 'lpv_c2', c2, 'lpv_c3', c3, 'lpv_c4', c4, ...
        'lpv_c5', c5, 'lpv_c6', c6, 'lpv_c7', c7, 'lpv_invV', invV, ...
        'lpv_Q', Q, 'lpv_V', V, 'lpv_h', h, ...
        'Kp_th0', K0.Kp_th, 'Kd_th0', K0.Kd_th, 'Kp_z0', K0.Kp_z, 'Kd_z0', K0.Kd_z, ...
        'tsched', tsched, 'Kp_sched', Kp_sched, 'Kd_sched', Kd_sched, ...
        'sched', 0, 'drywind', drywind, 'GreensiteLPV', GreensiteLPV, ...
        'Tstart', 0, 'Tstop', o.Tstop);
    fn = fieldnames(base);
    for i = 1:numel(fn), assignin('base', fn{i}, base.(fn{i})); end
    fprintf('init_simulink_lpv: pushed %d variables to base (horizon 0-%g s, %d schedule points).\n', ...
            numel(fn), o.Tstop, numel(tsched));
end
end

% ------------------------------------------------------------------------
function wg = run_wind_generator(gdir, Tstop, drywind, GreensiteLPV)
%RUN_WIND_GENERATOR  Full-ascent wind from strong_wind.slx (never modified).
%   The generator is loaded in memory, its two component outputs (mean
%   profile v_wp and altitude-scheduled Dryden turbulence) are marked for
%   logging, the model is simulated over [0, Tstop] with fixed internal noise
%   seeds (reproducible), and closed WITHOUT saving. Returns the total wind
%   v_w = v_wp + turbulence on the union of the two (variable-step) log grids.
load_system(fullfile(gdir, 'strong_wind.slx'));
cleanup = onCleanup(@() close_system('strong_wind', 0));        % discard edits

ph    = get_param('strong_wind/Subsystem', 'PortHandles');
names = {'sw_vwp', 'sw_turb'};
for k = 1:2                          % logging lives on the source PORT
    set_param(ph.Outport(k), 'Name', names{k}, 'DataLogging', 'on');
end

in = Simulink.SimulationInput('strong_wind');
in = in.setVariable('drywind', drywind);
in = in.setVariable('GreensiteLPV', GreensiteLPV);
in = in.setModelParameter('StopTime', num2str(Tstop), ...
                          'SignalLogging', 'on', 'SignalLoggingName', 'logsout');
so = sim(in);

vwp  = so.logsout.getElement('sw_vwp').Values;
turb = so.logsout.getElement('sw_turb').Values;
[t1u, i1] = unique(vwp.Time);        % variable-step logs may repeat times
[t2u, i2] = unique(turb.Time);
d1 = squeeze(vwp.Data);  d2 = squeeze(turb.Data);
t  = unique([t1u; t2u]);
vw = interp1(t1u, d1(i1), t, 'linear', 'extrap') + ...
     interp1(t2u, d2(i2), t, 'linear', 'extrap');
wg = struct('t', t(:), 'vw', vw(:));
end
