function w = load_wind_profile(p, varargin)
%LOAD_WIND_PROFILE  Wind angle-of-attack disturbance for the gust simulation.
%
%   w = LOAD_WIND_PROFILE(p) returns a struct describing a wind disturbance
%   to drive the closed-loop time response (Task 1-3). The disturbance is
%   expressed as a wind-induced angle of attack
%
%       alpha_w(t) = v_w(t) / V
%
%   where v_w is the lateral wind velocity [m/s] and V is the LV relative
%   speed at max-qbar. The default is a deterministic "1-cosine" discrete
%   gust whose amplitude is taken from the reference dry-wind dispersion
%   (drywind.mat) at the current altitude, severe regime.
%
%   Name/value options:
%     'profile'  'gust' (default) | 'step' | 'doublet' | 'strongwind'
%     'severity' 'light'|'moderate'|'severe' (default 'severe')
%     'Vg'       gust peak velocity [m/s]    (default: drywind sigma @ Alt)
%     'Tg'       gust duration [s]           (default 3 s)
%     'Tend'     simulation horizon [s]      (default 12 s)
%     'dt'       sample time [s]             (default 0.005 s)
%     't0'       gust onset time [s]         (default 1 s)
%
%   'strongwind' runs the professor's wind generator
%   (General/hw3-v3/strong_wind.slx: mean profile + altitude-scheduled
%   Dryden turbulence, fixed noise seed) and windows the total wind
%   v_wp + noisy_wind around the max-qbar instant p.t_ref, so the frozen
%   time closed loop is driven by the wind it would actually see there.
%   The .slx is loaded in memory and never modified.
%
%   Output fields: t (1xN), vw (1xN) [m/s], alphaw (1xN) [rad], plus the
%   metadata used to build it.
%
%   See also SIMULATE_GUST_RESPONSE, LOAD_HW3_PARAMS.

ip = inputParser;
ip.addParameter('profile','gust');
ip.addParameter('severity','severe');
ip.addParameter('Vg',[]);
ip.addParameter('Tg',3.0);
ip.addParameter('Tend',12.0);
ip.addParameter('dt',0.005);
ip.addParameter('t0',1.0);
ip.parse(varargin{:});
o = ip.Results;

% --- professor's generator: simulate strong_wind.slx and window at max-qbar ---
if strcmpi(o.profile, 'strongwind')
    w = local_strong_wind(p, o);
    return;
end

% --- default gust amplitude from drywind dispersion at current altitude ---
Vg = o.Vg;
if isempty(Vg)
    dwfile = fullfile(fileparts(mfilename('fullpath')), ...
                      'General','hw3-v3','drywind.mat');
    if isfile(dwfile)
        S = load(dwfile); dw = S.drywind;
        alt_km = p.Alt/1000;                       % drywind.alt is in km
        sig = dw.sigma.(o.severity);
        Vg = interp1(dw.alt, sig, alt_km, 'linear', 'extrap');
    else
        Vg = 8.0;                                  % fallback [m/s]
    end
end

t  = 0:o.dt:o.Tend;
vw = zeros(size(t));
switch lower(o.profile)
    case 'gust'   % 1-cosine discrete gust
        idx = t >= o.t0 & t <= o.t0 + o.Tg;
        vw(idx) = 0.5*Vg*(1 - cos(2*pi*(t(idx)-o.t0)/o.Tg));
    case 'step'
        vw(t >= o.t0) = Vg;
    case 'doublet'
        idx1 = t >= o.t0 & t < o.t0 + o.Tg/2;
        idx2 = t >= o.t0 + o.Tg/2 & t <= o.t0 + o.Tg;
        vw(idx1) =  Vg;  vw(idx2) = -Vg;
    otherwise
        error('load_wind_profile:profile','unknown profile ''%s''.',o.profile);
end

w = struct('t',t,'vw',vw,'alphaw',vw/p.V, 'V',p.V, ...
           'Vg',Vg,'Tg',o.Tg,'profile',o.profile,'severity',o.severity);
end

function w = local_strong_wind(p, o)
%LOCAL_STRONG_WIND  Total wind from strong_wind.slx, windowed at max-qbar.
%   The generator model is loaded in memory, its two outputs (mean profile
%   v_wp and Dryden turbulence) are marked for signal logging, and the model
%   is closed WITHOUT saving, so the professor's file is never touched. The
%   noise seeds are fixed inside the model, so the run is reproducible.
here = fileparts(mfilename('fullpath'));
gdir = fullfile(here,'General','hw3-v3');
dw   = load(fullfile(gdir,'drywind.mat'));
lpv  = load(fullfile(gdir,'GreensiteLPV_DATA.mat'));

load_system(fullfile(gdir,'strong_wind.slx'));
cleanup = onCleanup(@() close_system('strong_wind', 0));   % discard edits

ph = get_param('strong_wind/Subsystem','PortHandles');
names = {'sw_vwp','sw_turb'};
for k = 1:2                      % logging lives on the source PORT, not the line
    set_param(ph.Outport(k), 'Name', names{k}, 'DataLogging', 'on');
end

t1 = p.t_ref - o.t0;                       % window start, gust-onset aligned
in = Simulink.SimulationInput('strong_wind');
in = in.setVariable('drywind', dw.drywind);
in = in.setVariable('GreensiteLPV', lpv.GreensiteLPV);
in = in.setModelParameter('StopTime', num2str(t1 + o.Tend), ...
                          'SignalLogging','on', 'SignalLoggingName','logsout');
so = sim(in);

vwp  = so.logsout.getElement('sw_vwp').Values;
turb = so.logsout.getElement('sw_turb').Values;

t  = 0:o.dt:o.Tend;
[t1u, i1] = unique(vwp.Time);              % variable-step logs may repeat times
[t2u, i2] = unique(turb.Time);
d1 = squeeze(vwp.Data);  d2 = squeeze(turb.Data);
vw = interp1(t1u, d1(i1), t1 + t, 'linear') + ...
     interp1(t2u, d2(i2), t1 + t, 'linear');

w = struct('t',t,'vw',vw,'alphaw',vw/p.V, 'V',p.V, ...
           'Vg',max(abs(vw)),'Tg',o.Tend,'profile','strongwind', ...
           'severity','severe','t_window',[t1, t1 + o.Tend]);
end
