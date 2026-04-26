%% main.m - Falcon 9 First Stage Trajectory Simulation
%  Classwork n. 1 - Simulation of Falcon 9 First Stage Motion
%  3-DoF point-mass model in spherical coordinates
%  Velocity in Up-East-North (UEN) frame

clear; close all; clc;

%% ========================================================================
%  CONSTANTS AND VEHICLE PARAMETERS
%  ========================================================================

% --- Planetary constants ---
mu     = 3.986004418e14;     % Earth gravitational parameter [m^3/s^2]
RE     = 6378137;            % Earth equatorial radius [m]

% --- Atmospheric model (exponential, isothermal) ---
rho0   = 1.225;              % Sea-level air density [kg/m^3]
p0     = 101325;             % Sea-level pressure [Pa]
Hscale = 8000;               % Scale height [m] (8 km)
Tamb   = 288.15;             % Ambient temperature [K]

% --- Earth rotation ---
Tsid   = 86136;              % Sidereal day [s]
omegaE = 2 * pi / Tsid;     % Earth angular velocity [rad/s]

% --- Thermodynamic ---
gamma_air = 1.4;             % Ratio of specific heats [-]
Rgas      = 287.058;         % Specific gas constant for air [J/(kg*K)]
a_sound   = sqrt(gamma_air * Rgas * Tamb); % Speed of sound [m/s]

% --- First Stage ---
mdry1  = 22200;              % Dry mass [kg]
mp1    = 410900;             % Propellant mass [kg]
Tvac1  = 8227e3;             % Vacuum thrust [N]
tb1    = 162;                % Burn time [s]
cvac1  = 3244;               % Vacuum exhaust velocity [m/s]
Aex1   = 11.039;             % Nozzle exit area [m^2]
Qdot1  = Tvac1 / cvac1;     % Mass flow rate [kg/s]

% --- Second Stage (mass contribution only) ---
mdry2  = 4000;               % Dry mass [kg]
mp2    = 107500;             % Propellant mass [kg]

% --- Aerodynamics ---
CD     = 0.329;              % Drag coefficient [-]
Sref   = 10.52;              % Reference area [m^2]

% --- Payload and fairing ---
mfair  = 1700;               % Fairing mass [kg]
mpay   = 22800;              % Payload mass [kg] (LEO capacity)

% --- Phase timing ---
t_end1 = 5;                  % End of vertical ascent [s]
t_end2 = 15;                 % End of pitchover maneuver [s]

%% ========================================================================
%  INITIAL CONDITIONS - Launch from Kennedy Space Center
%  ========================================================================

lat0 = deg2rad(28.573469);   % Latitude [rad]
lon0 = deg2rad(-80.651070);  % Longitude [rad]

r0     = RE;                               % Radius [m]
theta0 = lon0;                             % Right ascension [rad]
phi0   = lat0;                             % Declination [rad]
u0     = 0;                                % Up velocity [m/s]
v0     = omegaE * RE * cos(lat0);          % East velocity (Earth rot.) [m/s]
w0     = 0;                                % North velocity [m/s]
m0     = mdry1 + mp1 + mdry2 + mp2 + mfair + mpay;  % Total mass [kg]

y0 = [r0; theta0; phi0; u0; v0; w0; m0];

%% ========================================================================
%  PARAMETER STRUCTURE
%  ========================================================================

par.mu     = mu;
par.RE     = RE;
par.rho0   = rho0;
par.p0     = p0;
par.H      = Hscale;
par.omegaE = omegaE;
par.Tvac   = Tvac1;
par.cvac   = cvac1;
par.Aex    = Aex1;
par.Qdot   = Qdot1;
par.CD     = CD;
par.Sref   = Sref;
par.t1     = t_end1;
par.t2     = t_end2;

%% ========================================================================
%  NUMERICAL INTEGRATION (ode45)
%  ========================================================================

opts  = odeset('RelTol', 1e-10, 'AbsTol', 1e-12, 'MaxStep', 1);
tspan = [0, tb1];
[t, Y] = ode45(@(t, y) eom(t, y, par), tspan, y0, opts);

%% ========================================================================
%  POST-PROCESSING
%  ========================================================================

% --- Extract state variables ---
r     = Y(:,1);
theta = Y(:,2);
phi   = Y(:,3);
uvel  = Y(:,4);
vvel  = Y(:,5);
wvel  = Y(:,6);
mass  = Y(:,7);

% --- Altitude ---
h = r - RE;

% --- Inertial velocity magnitude ---
Vmag = sqrt(uvel.^2 + vvel.^2 + wvel.^2);

% --- Relative velocity (atmosphere co-rotates with Earth) ---
urel = uvel;
vrel = vvel - omegaE .* r .* cos(phi);
wrel = wvel;
Vrel = sqrt(urel.^2 + vrel.^2 + wrel.^2);

% --- Atmospheric quantities ---
rho_traj = rho0 * exp(-h / Hscale);
p_traj   = p0   * exp(-h / Hscale);

% --- Dynamic pressure ---
qdyn = 0.5 * rho_traj .* Vrel.^2;

% --- Mach number ---
Mach = Vrel / a_sound;

% --- Thrust magnitude ---
Tmag = Tvac1 - p_traj * Aex1;

% --- Thrust elevation and aerodynamic flight-path angle ---
gammaT = zeros(size(t));
gammaA = zeros(size(t));

for k = 1:length(t)
    % Aerodynamic angle: elevation of relative velocity
    Vh = sqrt(vrel(k)^2 + wrel(k)^2);
    if Vrel(k) > 1e-6
        gammaA(k) = atan2d(urel(k), Vh);
    else
        gammaA(k) = 90;
    end

    % Thrust elevation angle
    if t(k) <= t_end1
        gammaT(k) = 90;                                  % Vertical ascent
    elseif t(k) <= t_end2
        gammaT(k) = 90 - 0.05 * (t(k) - t_end1);        % Pitchover
    else
        gammaT(k) = gammaA(k);                            % Gravity turn
    end
end

% --- Identify mission events ---
[~, ip1] = min(abs(t - t_end1));          % End of vertical ascent
[~, ip2] = min(abs(t - t_end2));          % End of pitchover
im1 = find(Mach >= 1, 1, 'first');        % Mach 1
[qmax, imQ] = max(qdyn);                  % Max-Q

% --- Local ENU displacements for 3D trajectory ---
% Cartesian ECI coordinates
Xeci = r .* cos(phi) .* cos(theta);
Yeci = r .* cos(phi) .* sin(theta);
Zeci = r .* sin(phi);

% Launch point in ECI
X0 = RE * cos(phi0) * cos(theta0);
Y0 = RE * cos(phi0) * sin(theta0);
Z0 = RE * sin(phi0);

% Rotation matrix ECI -> ENU at launch site
R_enu = [-sin(theta0),             cos(theta0),              0;
         -sin(phi0)*cos(theta0),  -sin(phi0)*sin(theta0),    cos(phi0);
          cos(phi0)*cos(theta0),   cos(phi0)*sin(theta0),    sin(phi0)];

dR_eci = [Xeci - X0, Yeci - Y0, Zeci - Z0];
dR_enu = (R_enu * dR_eci')';

East_km  = dR_enu(:,1) / 1000;
North_km = dR_enu(:,2) / 1000;
Up_km    = dR_enu(:,3) / 1000;

% --- Ground track (ECEF longitude) ---
lon_ground = rad2deg(theta - omegaE * t);
lat_ground = rad2deg(phi);

%% ========================================================================
%  CONSOLE SUMMARY
%  ========================================================================

fprintf('====================================================\n');
fprintf('  FALCON 9 FIRST STAGE - SIMULATION RESULTS\n');
fprintf('====================================================\n');
fprintf('  Initial mass:        %9.0f kg\n', m0);
fprintf('  Final mass:          %9.0f kg\n', mass(end));
fprintf('  Propellant consumed: %9.0f kg\n', m0 - mass(end));
fprintf('  Mass flow rate:      %9.2f kg/s\n', Qdot1);
fprintf('----------------------------------------------------\n');
fprintf('  Final altitude:      %9.2f km\n', h(end)/1000);
fprintf('  Final inertial vel.: %9.1f m/s\n', Vmag(end));
fprintf('  Final relative vel.: %9.1f m/s\n', Vrel(end));
fprintf('  Final Mach:          %9.2f\n', Mach(end));
fprintf('----------------------------------------------------\n');
if ~isempty(im1)
    fprintf('  Mach 1:    t = %6.1f s   h = %6.2f km\n', t(im1), h(im1)/1e3);
end
fprintf('  Max-Q:     t = %6.1f s   h = %6.2f km   q = %.1f kPa\n', ...
    t(imQ), h(imQ)/1e3, qmax/1e3);
fprintf('====================================================\n');

%% ========================================================================
%  PLOTS
%  ========================================================================

% Event marker colors
cP1 = [0.2 0.7 0.2];    % End vertical ascent (green)
cP2 = [0.9 0.6 0.0];    % End pitchover (orange)
cM1 = [0.85 0.0 0.0];   % Mach 1 (red)
cMQ = [0.0 0.2 0.8];    % Max-Q (blue)

% ---- Figure 1: Three-Dimensional Trajectory (Local ENU) ----
figure('Name', '3D Trajectory', 'Position', [50 400 700 550]);
plot3(East_km, North_km, Up_km, 'b', 'LineWidth', 1.8);
hold on;
plot3(East_km(1), North_km(1), Up_km(1), ...
    'go', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
plot3(East_km(end), North_km(end), Up_km(end), ...
    'rs', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
plot3(East_km(ip1), North_km(ip1), Up_km(ip1), ...
    'v', 'Color', cP1, 'MarkerSize', 10, 'MarkerFaceColor', cP1);
plot3(East_km(ip2), North_km(ip2), Up_km(ip2), ...
    'v', 'Color', cP2, 'MarkerSize', 10, 'MarkerFaceColor', cP2);
if ~isempty(im1)
    plot3(East_km(im1), North_km(im1), Up_km(im1), ...
        '^', 'Color', cM1, 'MarkerSize', 10, 'MarkerFaceColor', cM1);
end
plot3(East_km(imQ), North_km(imQ), Up_km(imQ), ...
    'd', 'Color', cMQ, 'MarkerSize', 10, 'MarkerFaceColor', cMQ);
xlabel('East [km]'); ylabel('North [km]'); zlabel('Up [km]');
title('Falcon 9 First Stage - 3D Trajectory');
legend('Trajectory', 'Launch', 'MECO', 'End Vert. Ascent', ...
    'End Pitchover', 'Mach 1', 'Max-Q', 'Location', 'best');
grid on; view([-35 25]); hold off;

% ---- Figure 2: Altitude vs Time ----
figure('Name', 'Altitude', 'Position', [100 350 700 450]);
plot(t, h/1e3, 'b', 'LineWidth', 1.5); hold on;
xline(t_end1, '--', 'Color', cP1, 'LineWidth', 1, 'HandleVisibility', 'off');
xline(t_end2, '--', 'Color', cP2, 'LineWidth', 1, 'HandleVisibility', 'off');
plot(t(ip1), h(ip1)/1e3, 'v', 'Color', cP1, 'MarkerSize', 9, ...
    'MarkerFaceColor', cP1, 'DisplayName', 'End Vert. Ascent');
plot(t(ip2), h(ip2)/1e3, 'v', 'Color', cP2, 'MarkerSize', 9, ...
    'MarkerFaceColor', cP2, 'DisplayName', 'End Pitchover');
if ~isempty(im1)
    plot(t(im1), h(im1)/1e3, '^', 'Color', cM1, 'MarkerSize', 9, ...
        'MarkerFaceColor', cM1, 'DisplayName', sprintf('Mach 1 (t=%.1fs)', t(im1)));
end
plot(t(imQ), h(imQ)/1e3, 'd', 'Color', cMQ, 'MarkerSize', 9, ...
    'MarkerFaceColor', cMQ, 'DisplayName', sprintf('Max-Q (t=%.1fs)', t(imQ)));
xlabel('Time [s]'); ylabel('Altitude [km]');
title('Altitude vs Time');
legend('Altitude', 'Location', 'northwest'); grid on; hold off;

% ---- Figure 3: Velocity Magnitude vs Time ----
figure('Name', 'Velocity', 'Position', [150 300 700 450]);
plot(t, Vmag, 'b', 'LineWidth', 1.5, 'DisplayName', 'Inertial |V|'); hold on;
plot(t, Vrel, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Relative |V_{rel}|');
xline(t_end1, '--', 'Color', cP1, 'LineWidth', 1, 'HandleVisibility', 'off');
xline(t_end2, '--', 'Color', cP2, 'LineWidth', 1, 'HandleVisibility', 'off');
if ~isempty(im1)
    plot(t(im1), Vrel(im1), '^', 'Color', cM1, 'MarkerSize', 9, ...
        'MarkerFaceColor', cM1, 'DisplayName', sprintf('Mach 1 (t=%.1fs)', t(im1)));
end
plot(t(imQ), Vrel(imQ), 'd', 'Color', cMQ, 'MarkerSize', 9, ...
    'MarkerFaceColor', cMQ, 'DisplayName', sprintf('Max-Q (t=%.1fs)', t(imQ)));
xlabel('Time [s]'); ylabel('Velocity [m/s]');
title('Velocity Magnitude vs Time');
legend('Location', 'northwest'); grid on; hold off;

% ---- Figure 4: Mass vs Time ----
figure('Name', 'Mass', 'Position', [200 250 700 450]);
plot(t, mass/1e3, 'b', 'LineWidth', 1.5); hold on;
xline(t_end1, '--', 'Color', cP1, 'LineWidth', 1, 'HandleVisibility', 'off');
xline(t_end2, '--', 'Color', cP2, 'LineWidth', 1, 'HandleVisibility', 'off');
plot(t(ip1), mass(ip1)/1e3, 'v', 'Color', cP1, 'MarkerSize', 9, ...
    'MarkerFaceColor', cP1, 'DisplayName', 'End Vert. Ascent');
plot(t(ip2), mass(ip2)/1e3, 'v', 'Color', cP2, 'MarkerSize', 9, ...
    'MarkerFaceColor', cP2, 'DisplayName', 'End Pitchover');
xlabel('Time [s]'); ylabel('Mass [t]');
title('Vehicle Mass vs Time');
legend('Mass', 'Location', 'northeast'); grid on; hold off;

% ---- Figure 5: Dynamic Pressure vs Time ----
figure('Name', 'Dynamic Pressure', 'Position', [250 200 700 450]);
plot(t, qdyn/1e3, 'b', 'LineWidth', 1.5, 'DisplayName', 'q_{dyn}'); hold on;
xline(t_end1, '--', 'Color', cP1, 'LineWidth', 1, 'HandleVisibility', 'off');
xline(t_end2, '--', 'Color', cP2, 'LineWidth', 1, 'HandleVisibility', 'off');
if ~isempty(im1)
    plot(t(im1), qdyn(im1)/1e3, '^', 'Color', cM1, 'MarkerSize', 9, ...
        'MarkerFaceColor', cM1, 'DisplayName', sprintf('Mach 1 (t=%.1fs)', t(im1)));
end
plot(t(imQ), qmax/1e3, 'd', 'Color', cMQ, 'MarkerSize', 12, ...
    'MarkerFaceColor', cMQ, 'DisplayName', sprintf('Max-Q = %.1f kPa', qmax/1e3));
xlabel('Time [s]'); ylabel('Dynamic Pressure [kPa]');
title('Dynamic Pressure vs Time');
legend('Location', 'northeast'); grid on; hold off;

% ---- Figure 6: Mach Number vs Time ----
figure('Name', 'Mach Number', 'Position', [300 150 700 450]);
plot(t, Mach, 'b', 'LineWidth', 1.5, 'DisplayName', 'Mach'); hold on;
yline(1, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
xline(t_end1, '--', 'Color', cP1, 'LineWidth', 1, 'HandleVisibility', 'off');
xline(t_end2, '--', 'Color', cP2, 'LineWidth', 1, 'HandleVisibility', 'off');
if ~isempty(im1)
    plot(t(im1), 1, '^', 'Color', cM1, 'MarkerSize', 9, ...
        'MarkerFaceColor', cM1, 'DisplayName', sprintf('Mach 1 (t=%.1fs)', t(im1)));
end
xlabel('Time [s]'); ylabel('Mach Number [-]');
title('Mach Number vs Time');
legend('Location', 'northwest'); grid on; hold off;

% ---- Figure 7: Thrust and Aerodynamic Angles vs Time ----
figure('Name', 'Angles', 'Position', [350 100 700 450]);
plot(t, gammaT, 'b', 'LineWidth', 1.5, 'DisplayName', '\gamma_T (Thrust)'); hold on;
plot(t, gammaA, 'r--', 'LineWidth', 1.5, 'DisplayName', '\gamma_A (Aerodynamic)');
xline(t_end1, '--', 'Color', cP1, 'LineWidth', 1, 'HandleVisibility', 'off');
xline(t_end2, '--', 'Color', cP2, 'LineWidth', 1, 'HandleVisibility', 'off');
xlabel('Time [s]'); ylabel('Angle [deg]');
title('Thrust Elevation and Aerodynamic Flight-Path Angle vs Time');
legend('Location', 'best'); grid on; hold off;

% ---- Figure 8: Ground Track ----
figure('Name', 'Ground Track', 'Position', [400 50 700 450]);
plot(lon_ground, lat_ground, 'b', 'LineWidth', 2); hold on;
plot(lon_ground(1), lat_ground(1), 'go', 'MarkerSize', 12, ...
    'MarkerFaceColor', 'g', 'DisplayName', 'Launch (KSC)');
plot(lon_ground(end), lat_ground(end), 'rs', 'MarkerSize', 12, ...
    'MarkerFaceColor', 'r', 'DisplayName', 'MECO');
if ~isempty(im1)
    plot(lon_ground(im1), lat_ground(im1), '^', 'Color', cM1, ...
        'MarkerSize', 9, 'MarkerFaceColor', cM1, 'DisplayName', 'Mach 1');
end
plot(lon_ground(imQ), lat_ground(imQ), 'd', 'Color', cMQ, ...
    'MarkerSize', 9, 'MarkerFaceColor', cMQ, 'DisplayName', 'Max-Q');
xlabel('Longitude [deg]'); ylabel('Latitude [deg]');
title('Ground Track');
legend('Trajectory', 'Location', 'best'); grid on; hold off;

%% ========================================================================
%  EXPORT FIGURES TO figures/ (PNG, 200 dpi)
%  ========================================================================
fig_dir = fullfile(fileparts(mfilename('fullpath')), 'figures');
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end

slugify = @(s) lower(regexprep(s, '[^a-zA-Z0-9]+', '_'));
fig_handles = findobj(groot, 'Type', 'figure');
for kk = 1:numel(fig_handles)
    fname = fullfile(fig_dir, [slugify(get(fig_handles(kk), 'Name')) '.png']);
    exportgraphics(fig_handles(kk), fname, 'Resolution', 200);
end

%% ========================================================================
%  EQUATIONS OF MOTION (local function)
%  ========================================================================

function dydt = eom(t, y, par)
% 3-DoF point-mass equations in spherical coordinates (r, theta, phi)
% with inertial velocity in Up-East-North frame: V = [u, v, w]
%
% State vector: y = [r, theta, phi, u, v, w, m]
%
% Forces: Keplerian gravity, aerodynamic drag, rocket thrust
% Three flight phases:
%   1) Vertical ascent        (0        <= t <= t1)
%   2) Pitchover maneuver     (t1       <  t <= t2)
%   3) Zero-lift gravity turn (t2       <  t <= tb)

    % --- Unpack state ---
    r     = y(1);   % Radius [m]
    phi   = y(3);   % Declination / latitude [rad]
    u     = y(4);   % Up velocity [m/s]
    v     = y(5);   % East velocity [m/s]
    w     = y(6);   % North velocity [m/s]
    m     = y(7);   % Mass [kg]

    % --- Altitude ---
    alt = r - par.RE;

    % --- Atmosphere (exponential model) ---
    rho  = par.rho0 * exp(-alt / par.H);
    patm = par.p0   * exp(-alt / par.H);

    % --- Relative velocity (atmosphere co-rotates) ---
    v_atm = par.omegaE * r * cos(phi);
    ur = u;
    vr = v - v_atm;
    wr = w;
    Vrel = sqrt(ur^2 + vr^2 + wr^2);

    % --- Gravity (inverse-square, radial only) ---
    g_u = -par.mu / r^2;

    % --- Aerodynamic drag (opposite to Vrel) ---
    if Vrel > 1e-10
        Dcoeff = -0.5 * rho * par.Sref * par.CD * Vrel;
        Du = Dcoeff * ur;
        Dv = Dcoeff * vr;
        Dw = Dcoeff * wr;
    else
        Du = 0; Dv = 0; Dw = 0;
    end

    % --- Thrust ---
    Tmag = par.Tvac - patm * par.Aex;

    if t <= par.t1
        % Phase 1: Vertical ascent - thrust purely upward
        Tu = Tmag;
        Tv = 0;
        Tw = 0;

    elseif t <= par.t2
        % Phase 2: Pitchover maneuver
        %   Elevation: linear 90 deg -> 89.5 deg over 10 s
        %   Azimuth:   constant 90 deg (East)
        gT  = deg2rad(90 - 0.05 * (t - par.t1));
        psi = deg2rad(90);
        Tu = Tmag * sin(gT);
        Tv = Tmag * cos(gT) * sin(psi);
        Tw = Tmag * cos(gT) * cos(psi);

    else
        % Phase 3: Zero-lift gravity turn - thrust along Vrel
        if Vrel > 1e-10
            Tu = Tmag * ur / Vrel;
            Tv = Tmag * vr / Vrel;
            Tw = Tmag * wr / Vrel;
        else
            Tu = Tmag; Tv = 0; Tw = 0;
        end
    end

    % --- Kinematic equations ---
    dr     = u;
    dtheta = v / (r * cos(phi));
    dphi   = w / r;

    % --- Dynamic equations (transport terms + forces/m) ---
    du = (v^2 + w^2) / r  +  g_u  +  (Tu + Du) / m;
    dv = (-u*v + v*w*tan(phi)) / r  +  (Tv + Dv) / m;
    dw = (-u*w - v^2*tan(phi)) / r  +  (Tw + Dw) / m;

    % --- Mass equation ---
    dm = -par.Qdot;

    dydt = [dr; dtheta; dphi; du; dv; dw; dm];
end