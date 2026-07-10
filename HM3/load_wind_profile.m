function w = load_wind_profile(p, o)
% Wind angle-of-attack disturbance for the gust simulation. alpha_w = v_w/V.
%   INPUT
%     p - param struct (uses V, Alt, t_ref)
%     o - name-value:
%           profile  'gust'(default)|'step'|'doublet'|'strongwind'
%           severity 'light'|'moderate'|'severe' (default 'severe')
%           Vg       gust peak velocity [m/s]    (default: drywind sigma @ Alt)
%           Tg       gust duration [s]           (default 3)
%           Tend     horizon [s]                 (default 12)
%           dt       sample time [s]             (default 0.005)
%           t0       gust onset [s]              (default 1)
%   OUTPUT
%     w - struct: t, vw [m/s], alphaw [rad], plus build metadata
%   'gust' = 1-cosine discrete gust, amplitude from drywind.mat dispersion.
%   'strongwind' runs strong_wind.slx (mean + altitude-scheduled Dryden, fixed
%   seed) and windows the total wind at max-qbar; the .slx is never modified.

arguments
    p (1,1) struct
    o.profile  {mustBeTextScalar} = 'gust'   % case-insensitive below
    o.severity {mustBeMember(o.severity, ["light","moderate","severe"])} = 'severe'
    o.Vg   {mustBeNumeric, mustBeReal, mustBeScalarOrEmpty} = []
    o.Tg   (1,1) {mustBeNumeric, mustBeReal, mustBePositive} = 3.0
    o.Tend (1,1) {mustBeNumeric, mustBeReal, mustBePositive} = 12.0
    o.dt   (1,1) {mustBeNumeric, mustBeReal, mustBePositive} = 0.005
    o.t0   (1,1) {mustBeNumeric, mustBeReal, mustBeNonnegative} = 1.0
end

% --- strongwind: simulate strong_wind.slx and window at max-qbar ---
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
% Total wind from strong_wind.slx, windowed at max-qbar.
%   INPUT  p, o - as in load_wind_profile
%   OUTPUT w    - wind struct (same shape as the analytic profiles)
%   Model loaded in memory, its two outputs (mean v_wp, Dryden turbulence)
%   logged, closed WITHOUT saving (file untouched). Seeds fixed -> reproducible.
here = fileparts(mfilename('fullpath'));
gdir = fullfile(here,'General','hw3-v3');
dw   = load(fullfile(gdir,'drywind.mat'));
lpv  = load(fullfile(gdir,'GreensiteLPV_DATA.mat'));

load_system(fullfile(gdir,'strong_wind.slx'));
cleanup = onCleanup(@() close_system('strong_wind', 0));   % discard edits

ph = get_param('strong_wind/Subsystem','PortHandles');
names = {'sw_vwp','sw_turb'};
for k = 1:2                      % logging is on the source PORT, not the line
    set_param(ph.Outport(k), 'Name', names{k}, 'DataLogging', 'on');
end

t1 = p.t_ref - o.t0;                       % window start, aligned to gust onset
in = Simulink.SimulationInput('strong_wind');
in = in.setVariable('drywind', dw.drywind);
in = in.setVariable('GreensiteLPV', lpv.GreensiteLPV);
in = in.setModelParameter('StopTime', num2str(t1 + o.Tend), ...
                          'SignalLogging','on', 'SignalLoggingName','logsout');
so = sim(in);

vwp  = so.logsout.getElement('sw_vwp').Values;
turb = so.logsout.getElement('sw_turb').Values;

t  = 0:o.dt:o.Tend;
[t1u, i1] = unique(vwp.Time);              % variable-step logs can repeat times
[t2u, i2] = unique(turb.Time);
d1 = squeeze(vwp.Data);  d2 = squeeze(turb.Data);
vw = interp1(t1u, d1(i1), t1 + t, 'linear') + ...
     interp1(t2u, d2(i2), t1 + t, 'linear');

w = struct('t',t,'vw',vw,'alphaw',vw/p.V, 'V',p.V, ...
           'Vg',max(abs(vw)),'Tg',o.Tend,'profile','strongwind', ...
           'severity','severe','t_window',[t1, t1 + o.Tend]);
end
