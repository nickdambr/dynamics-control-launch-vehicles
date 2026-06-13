function p = load_hw3_params(opt)
%LOAD_HW3_PARAMS  Launch-vehicle parameters at the max-qbar condition (t = 72 s).
%
%   p = LOAD_HW3_PARAMS() returns a struct with the pitch-plane model
%   parameters of the Greensite fictitious launch vehicle at t = 72 s,
%   following Table 1 of "Homework 3 - Attitude Control of a Launch Vehicle
%   in Atmospheric Flight" (Zavoli, v1.2, 18 May 2026).
%
%   The time-varying coefficients are read from the reference data set
%   GreensiteLPV_DATA.mat (interpolated at t = 72 s) when available; the
%   geometry/mass quantities that only appear in Table 1 are hard-coded.
%   If the data file is missing, the whole table falls back to literals.
%
%   p = LOAD_HW3_PARAMS('mu_alpha_scale', sa, 'mu_c_scale', sc) scales the
%   aerodynamic moment coefficient mu_alpha = A6 and the control
%   effectiveness mu_c = K1 by the given factors (used in Task 3 for the
%   +/-30 % corner cases). Defaults are 1.0.
%
%   Field summary (units in comments):
%     A6,K1,a1,a3,a4   pitch/lateral coefficients   [1/s^2]
%     V                relative velocity            [m/s]
%     wBM,zBM          first bending mode freq/damp [rad/s], [-]
%     phi_ins,sigma_ins INS bending observation     [-], [rad/m]
%     phi_tvc          TVC bending forcing          [1/kg]
%     wTVC,zTVC,tau    TVC actuator + pure delay    [rad/s],[-],[s]

%% Options (name-value)
arguments
    opt.mu_alpha_scale (1,1) {mustBeNumeric, mustBeReal} = 1.0
    opt.mu_c_scale     (1,1) {mustBeNumeric, mustBeReal} = 1.0
    opt.t_ref          (1,1) {mustBeNumeric, mustBeReal} = 72
end

%% Table 1 literals (authoritative source: the assignment PDF)
p = struct();
p.t_ref         = opt.t_ref;     % s    reference instant (max-qbar)
% --- geometry / mass / forces (Table 1, time-invariant in this snapshot) ---
p.m             = 7.38e4;        % kg
p.l_alpha       = 10.39;         % m
p.l_c           = 9.84;          % m
p.Iyy           = 3.28e6;        % kg m^2
p.Alt           = 15143;         % m
p.Tt_minus_D    = 1.71e6;        % N
p.N_alpha       = 1.07e6;        % N/rad
% NOTE: Table 1 is internally inconsistent on a4: -(Tt-D)/m = -23.17 1/s^2
% with the Tt-D and m above, yet the table lists a4 = -27.2710. The LPV data
% set agrees with the latter, so a4 (not Tt-D) is what enters the dynamics.
% --- aero/control coefficients (Table 1 nominal) ---
p.A6            = 3.3818;        % 1/s^2  aerodynamic moment (mu_alpha)
p.K1            = 4.5647;        % 1/s^2  control effectiveness (mu_c)
p.a1            = -0.0154;       % 1/s^2
p.a3            = 20.6090;       % 1/s^2
p.a4            = -27.2710;      % 1/s^2
p.V             = 937.70;        % m/s
p.Tc            = 1.52e6;        % N
% --- bending mode ---
p.wBM           = 18.9;          % rad/s
p.zBM           = 0.005;         % -
% --- sensor (INS) model, Eq. (2) ---
p.phi_ins       = 0.8;           % -
p.sigma_ins     = 0.178;         % rad/m
% --- TVC actuator, Eq. (3) ---
p.phi_tvc       = 4.31e-5;       % 1/kg   bending forcing per unit thrust
p.wTVC          = 70;            % rad/s
p.zTVC          = 0.7;           % -
p.tau           = 0.020;         % s      pure transport delay

%% Prefer the reference LPV data set for the time-varying coefficients
datafile = fullfile(fileparts(mfilename('fullpath')), ...
                    'General', 'hw3-v3', 'GreensiteLPV_DATA.mat');
if isfile(datafile)
    S = load(datafile);
    L = S.GreensiteLPV;
    at = @(ts) interp1(ts.Time, squeeze(ts.Data), opt.t_ref);
    p.A6        = at(L.A6);
    p.K1        = at(L.K1);
    p.a1        = at(L.a1);
    p.a3        = at(L.a3);
    p.a4        = at(L.a4);
    p.V         = at(L.V);
    p.Tc        = at(L.Tc);
    p.sigma_ins = at(L.sigma_ins);
    p.phi_ins   = at(L.phi_ins);
    p.phi_tvc   = at(L.phi_tvc);
    p.wBM       = at(L.omega);
    p.src = 'GreensiteLPV_DATA.mat @ t=72 s';
else
    p.src = 'Table 1 literals (data file not found)';
end

%% Derived: dynamic pressure at the reference altitude (exponential atmosphere)
p.rho  = 1.225 * exp(-p.Alt/8000);  % kg/m^3
p.qbar = 0.5 * p.rho * p.V^2;       % Pa     (~81 kPa at max-qbar)

%% Apply Task-3 uncertainty scaling on the nominal coefficients
p.mu_alpha_scale = opt.mu_alpha_scale;
p.mu_c_scale     = opt.mu_c_scale;
p.A6 = p.A6 * opt.mu_alpha_scale;   % mu_alpha
p.K1 = p.K1 * opt.mu_c_scale;       % mu_c

end
