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
%     'profile'  'gust' (default) | 'step' | 'doublet'
%     'severity' 'light'|'moderate'|'severe' (default 'severe')
%     'Vg'       gust peak velocity [m/s]    (default: drywind sigma @ Alt)
%     'Tg'       gust duration [s]           (default 3 s)
%     'Tend'     simulation horizon [s]      (default 12 s)
%     'dt'       sample time [s]             (default 0.005 s)
%     't0'       gust onset time [s]         (default 1 s)
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
