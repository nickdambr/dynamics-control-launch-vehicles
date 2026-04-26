%% main2.m - Falcon 9 First Stage Trajectory Simulation (Non-Dimensional)
%  Classwork n. 1 - Simulation of Falcon 9 First Stage Motion
%  3-DoF point-mass model in spherical coordinates
%  Velocity in Up-East-North (UEN) frame
%
%  Non-dimensionalization:
%    L_ref = RE  = 6378137 m      (Earth equatorial radius)
%    V_ref       = 7800 m/s       (first cosmic velocity, ~7.8 km/s)
%    T_ref = L_ref / V_ref        (time reference)
%    m_ref = m0                   (initial total mass)
%
%  Derived reference scales:
%    a_ref = V_ref^2 / L_ref      (acceleration)
%    F_ref = m_ref * V_ref^2 / L_ref  (force)
%    p_ref = m_ref * V_ref^2 / L_ref^3 (pressure)
%
%  Integration: each flight arc is integrated over tau in [0, 1], where
%    tau = (t_dim - t_arc_start) / arc_duration
%  The EOM RHS is scaled by the non-dimensional arc duration Delta_k so
%  that  dy_nd/dtau = Delta_k * (dy_nd/dt_nd).

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
Hscale = 8000;               % Scale height [m]
Tamb   = 288.15;             % Ambient temperature [K]

% --- Earth rotation ---
Tsid   = 86136;              % Sidereal day [s]
omegaE = 2 * pi / Tsid;     % Earth angular velocity [rad/s]

% --- Thermodynamic ---
gamma_air = 1.4;
Rgas      = 287.058;
a_sound   = sqrt(gamma_air * Rgas * Tamb);  % Speed of sound [m/s]

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
mpay   = 22800;              % Payload mass [kg]

% --- Phase timing ---
t_end1 = 5;                  % End of vertical ascent [s]
t_end2 = 15;                 % End of pitchover maneuver [s]

%% ========================================================================
%  LAUNCH SITE - Kennedy Space Center
%  ========================================================================

lat0 = deg2rad(28.573469);   % Latitude  [rad]
lon0 = deg2rad(-80.651070);  % Longitude [rad]

%% ========================================================================
%  REFERENCE QUANTITIES AND NON-DIMENSIONALIZATION
%  ========================================================================

% Reference scales
L_ref = RE;                    % [m]
V_ref = 7800;                  % [m/s]  (first cosmic velocity ~7.8 km/s)
T_ref = L_ref / V_ref;         % [s]

% Mass reference = initial total mass (computed here before ICs)
m0    = mdry1 + mp1 + mdry2 + mp2 + mfair + mpay;
m_ref = m0;                    % [kg]

fprintf('=== Reference Scales ===\n');
fprintf('  L_ref = %.4f  km\n', L_ref/1e3);
fprintf('  V_ref = %.1f  m/s\n', V_ref);
fprintf('  T_ref = %.4f  s\n',  T_ref);
fprintf('  m_ref = %.0f  kg\n', m_ref);
fprintf('========================\n\n');

% --- Non-dimensional parameters ---
%   mu_nd    = mu  * T_ref^2 / L_ref^3  = mu / (V_ref^2 * L_ref)
%   omega_nd = omegaE * T_ref
%   rho0_nd  = rho0 * L_ref^3 / m_ref
%   p0_nd    = p0   * L_ref^3 / (m_ref * V_ref^2)  = p0 / p_ref
%   H_nd     = Hscale / L_ref
%   Tvac_nd  = Tvac1 * L_ref / (m_ref * V_ref^2)   = Tvac1 / F_ref
%   Aex_nd   = Aex1  / L_ref^2
%   Qdot_nd  = Qdot1 * T_ref / m_ref
%   Sref_nd  = Sref  / L_ref^2
%   CD       unchanged (dimensionless)

mu_nd      = mu     / (V_ref^2 * L_ref);
omegaE_nd  = omegaE * T_ref;
rho0_nd    = rho0   * L_ref^3 / m_ref;
p0_nd      = p0     * L_ref^3 / (m_ref * V_ref^2);
H_nd       = Hscale / L_ref;
Tvac_nd    = Tvac1  * L_ref  / (m_ref * V_ref^2);
Aex_nd     = Aex1   / L_ref^2;
Qdot_nd    = Qdot1  * T_ref  / m_ref;
Sref_nd    = Sref   / L_ref^2;

%% ========================================================================
%  NON-DIMENSIONAL INITIAL CONDITIONS
%  ========================================================================

r0_nd = RE / L_ref;             % = 1  (launch from Earth surface)
u0_nd = 0;
v0_nd = omegaE * RE * cos(lat0) / V_ref;  % Earth-rotation East velocity
w0_nd = 0;
m0_nd = 1;                      % m0 / m_ref

y0_nd = [r0_nd; lon0; lat0; u0_nd; v0_nd; w0_nd; m0_nd];

%% ========================================================================
%  ARC DURATIONS IN NON-DIMENSIONAL TIME
%  ========================================================================
%  Global nd time: t_nd = t_dim / T_ref
%  Arc 1: t_dim in [0,      t_end1]  -> t_nd in [0,    t1_nd]
%  Arc 2: t_dim in [t_end1, t_end2]  -> t_nd in [t1_nd, t2_nd]
%  Arc 3: t_dim in [t_end2, tb1  ]   -> t_nd in [t2_nd, tb_nd]

t1_nd  = t_end1 / T_ref;
t2_nd  = t_end2 / T_ref;
tb_nd  = tb1    / T_ref;

Delta1 = t1_nd;           % nd duration of arc 1
Delta2 = t2_nd - t1_nd;   % nd duration of arc 2
Delta3 = tb_nd - t2_nd;   % nd duration of arc 3

%% ========================================================================
%  PARAMETER STRUCTURE (all non-dimensional where applicable)
%  ========================================================================

par.mu       = mu_nd;
par.omegaE   = omegaE_nd;
par.rho0     = rho0_nd;
par.p0       = p0_nd;
par.H        = H_nd;
par.Tvac     = Tvac_nd;
par.Aex      = Aex_nd;
par.Qdot     = Qdot_nd;
par.CD       = CD;
par.Sref     = Sref_nd;
par.T_ref    = T_ref;     % [s]  needed to recover dim. time inside pitchover
par.t_end1   = t_end1;   % [s]  pitchover reference (dimensional)

%% ========================================================================
%  NUMERICAL INTEGRATION - THREE ARCS, EACH OVER tau in [0, 1]
%  ========================================================================
%  For arc k: tau = (t_dim - t_start_k) / duration_k
%  EOM:  dy_nd/dtau = Delta_k * (dy_nd/dt_nd)

opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);

% --- Arc 1: Vertical Ascent ---
par.phase       = 1;
par.Delta       = Delta1;
par.t_arc_start = 0;           % global nd time at start of arc
[tau1, Y1_nd] = ode45(@(tau,y) eom_nd(tau,y,par), [0 1], y0_nd, opts);

% --- Arc 2: Pitchover ---
par.phase       = 2;
par.Delta       = Delta2;
par.t_arc_start = t1_nd;
[tau2, Y2_nd] = ode45(@(tau,y) eom_nd(tau,y,par), [0 1], Y1_nd(end,:)', opts);

% --- Arc 3: Gravity Turn ---
par.phase       = 3;
par.Delta       = Delta3;
par.t_arc_start = t2_nd;
[tau3, Y3_nd] = ode45(@(tau,y) eom_nd(tau,y,par), [0 1], Y2_nd(end,:)', opts);

%% ========================================================================
%  RECONSTRUCT GLOBAL TIME AXIS AND COMBINE ARCS
%  ========================================================================

% Global non-dimensional time for each arc
t_nd_1 = tau1           * Delta1;
t_nd_2 = t1_nd + tau2   * Delta2;
t_nd_3 = t2_nd + tau3   * Delta3;

% Concatenate (drop duplicate boundary points)
t_nd = [t_nd_1;         t_nd_2(2:end);      t_nd_3(2:end)];
Y_nd = [Y1_nd;          Y2_nd(2:end,:);     Y3_nd(2:end,:)];

% Dimensional time [s]
t = t_nd * T_ref;

%% ========================================================================
%  POST-PROCESSING  (convert back to dimensional quantities)
%  ========================================================================

% State variables
r     = Y_nd(:,1) * L_ref;    % Radius [m]
theta = Y_nd(:,2);             % Right ascension [rad]
phi   = Y_nd(:,3);             % Declination / latitude [rad]
uvel  = Y_nd(:,4) * V_ref;    % Up velocity [m/s]
vvel  = Y_nd(:,5) * V_ref;    % East velocity [m/s]
wvel  = Y_nd(:,6) * V_ref;    % North velocity [m/s]
mass  = Y_nd(:,7) * m_ref;    % Mass [kg]

% Altitude
h = r - RE;

% Inertial velocity magnitude
Vmag = sqrt(uvel.^2 + vvel.^2 + wvel.^2);

% Relative velocity (atmosphere co-rotates with Earth)
urel = uvel;
vrel = vvel - omegaE .* r .* cos(phi);
wrel = wvel;
Vrel = sqrt(urel.^2 + vrel.^2 + wrel.^2);

% Atmospheric quantities along trajectory
rho_traj = rho0 * exp(-h / Hscale);
p_traj   = p0   * exp(-h / Hscale);

% Dynamic pressure
qdyn = 0.5 * rho_traj .* Vrel.^2;

% Mach number
Mach = Vrel / a_sound;

% Thrust magnitude
Tmag_traj = Tvac1 - p_traj * Aex1;

% Thrust elevation and aerodynamic flight-path angle
gammaT = zeros(size(t));
gammaA = zeros(size(t));

for k = 1:length(t)
    Vh = sqrt(vrel(k)^2 + wrel(k)^2);
    if Vrel(k) > 1e-6
        gammaA(k) = atan2d(urel(k), Vh);
    else
        gammaA(k) = 90;
    end
    if t(k) <= t_end1
        gammaT(k) = 90;
    elseif t(k) <= t_end2
        gammaT(k) = 90 - 0.05 * (t(k) - t_end1);
    else
        gammaT(k) = gammaA(k);
    end
end

% --- Identify mission events ---
[~, ip1] = min(abs(t - t_end1));         % End of vertical ascent
[~, ip2] = min(abs(t - t_end2));         % End of pitchover
im1      = find(Mach >= 1, 1, 'first');  % Mach 1
[qmax, imQ] = max(qdyn);                 % Max-Q

% --- ECI coordinates -> local ENU displacements ---
Xeci = r .* cos(phi) .* cos(theta);
Yeci = r .* cos(phi) .* sin(theta);
Zeci = r .* sin(phi);

X0 = RE * cos(lat0) * cos(lon0);
Y0 = RE * cos(lat0) * sin(lon0);
Z0 = RE * sin(lat0);

R_enu = [-sin(lon0),             cos(lon0),              0;
         -sin(lat0)*cos(lon0),  -sin(lat0)*sin(lon0),    cos(lat0);
          cos(lat0)*cos(lon0),   cos(lat0)*sin(lon0),    sin(lat0)];

dR_eci = [Xeci - X0, Yeci - Y0, Zeci - Z0];
dR_enu = (R_enu * dR_eci')';

East_km  = dR_enu(:,1) / 1e3;
North_km = dR_enu(:,2) / 1e3;
Up_km    = dR_enu(:,3) / 1e3;

% --- Ground track (ECEF longitude) ---
lon_ground = rad2deg(theta - omegaE * t);
lat_ground = rad2deg(phi);

%% ========================================================================
%  CONSOLE SUMMARY
%  ========================================================================

fprintf('====================================================\n');
fprintf('  FALCON 9 FIRST STAGE - SIMULATION RESULTS (ND)\n');
fprintf('====================================================\n');
fprintf('  Initial mass:        %9.0f kg\n', m0);
fprintf('  Final mass:          %9.0f kg\n', mass(end));
fprintf('  Propellant consumed: %9.0f kg\n', m0 - mass(end));
fprintf('  Mass flow rate:      %9.2f kg/s\n', Qdot1);
fprintf('----------------------------------------------------\n');
fprintf('  Final altitude:      %9.2f km\n', h(end)/1e3);
fprintf('  Final inertial vel.: %9.1f m/s\n', Vmag(end));
fprintf('  Final relative vel.: %9.1f m/s\n', Vrel(end));
fprintf('  Final Mach:          %9.2f\n',    Mach(end));
fprintf('----------------------------------------------------\n');
if ~isempty(im1)
    fprintf('  Mach 1:    t = %6.1f s   h = %6.2f km\n', t(im1), h(im1)/1e3);
end
fprintf('  Max-Q:     t = %6.1f s   h = %6.2f km   q = %.1f kPa\n', ...
    t(imQ), h(imQ)/1e3, qmax/1e3);
fprintf('----------------------------------------------------\n');
fprintf('  Non-dimensional final state:\n');
fprintf('    r* = %.8f       (r / RE)\n',    Y_nd(end,1));
fprintf('    u* = %.8f       (u / V_ref)\n', Y_nd(end,4));
fprintf('    v* = %.8f       (v / V_ref)\n', Y_nd(end,5));
fprintf('    w* = %.8f       (w / V_ref)\n', Y_nd(end,6));
fprintf('    m* = %.8f       (m / m0)\n',    Y_nd(end,7));
fprintf('    t* = %.8f       (t / T_ref)\n', t_nd(end));
fprintf('====================================================\n');

%% ========================================================================
%  PLOTS  (axes in dimensional units for readability)
%  ========================================================================

cP1 = [0.2 0.7 0.2];   % End vertical ascent (green)
cP2 = [0.9 0.6 0.0];   % End pitchover (orange)
cM1 = [0.85 0.0 0.0];  % Mach 1 (red)
cMQ = [0.0 0.2 0.8];   % Max-Q (blue)

% ---- Figure 1: Three-Dimensional Trajectory (Local ENU) ----
figure('Name', '3D Trajectory [ND sim]', 'Position', [50 400 700 550]);
plot3(East_km, North_km, Up_km, 'b', 'LineWidth', 1.8); hold on;
plot3(East_km(1),   North_km(1),   Up_km(1),   'go', 'MarkerSize', 12, ...
    'MarkerFaceColor', 'g');
plot3(East_km(end), North_km(end), Up_km(end), 'rs', 'MarkerSize', 12, ...
    'MarkerFaceColor', 'r');
plot3(East_km(ip1), North_km(ip1), Up_km(ip1), 'v', 'Color', cP1, ...
    'MarkerSize', 10, 'MarkerFaceColor', cP1);
plot3(East_km(ip2), North_km(ip2), Up_km(ip2), 'v', 'Color', cP2, ...
    'MarkerSize', 10, 'MarkerFaceColor', cP2);
if ~isempty(im1)
    plot3(East_km(im1), North_km(im1), Up_km(im1), '^', 'Color', cM1, ...
        'MarkerSize', 10, 'MarkerFaceColor', cM1);
end
plot3(East_km(imQ), North_km(imQ), Up_km(imQ), 'd', 'Color', cMQ, ...
    'MarkerSize', 10, 'MarkerFaceColor', cMQ);
xlabel('East [km]'); ylabel('North [km]'); zlabel('Up [km]');
title('Falcon 9 First Stage - 3D Trajectory (ND simulation)');
legend('Trajectory','Launch','MECO','End Vert. Ascent','End Pitchover', ...
    'Mach 1','Max-Q','Location','best');
grid on; view([-35 25]); hold off;

% ---- Figure 2: Altitude vs Time ----
figure('Name', 'Altitude [ND sim]', 'Position', [100 350 700 450]);
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
title('Altitude vs Time  (ND simulation)');
legend('Altitude', 'Location', 'northwest'); grid on; hold off;

% ---- Figure 3: Velocity Magnitude vs Time ----
figure('Name', 'Velocity [ND sim]', 'Position', [150 300 700 450]);
plot(t, Vmag, 'b',  'LineWidth', 1.5, 'DisplayName', 'Inertial |V|'); hold on;
plot(t, Vrel, 'r--','LineWidth', 1.5, 'DisplayName', 'Relative |V_{rel}|');
xline(t_end1, '--', 'Color', cP1, 'LineWidth', 1, 'HandleVisibility', 'off');
xline(t_end2, '--', 'Color', cP2, 'LineWidth', 1, 'HandleVisibility', 'off');
if ~isempty(im1)
    plot(t(im1), Vrel(im1), '^', 'Color', cM1, 'MarkerSize', 9, ...
        'MarkerFaceColor', cM1, 'DisplayName', sprintf('Mach 1 (t=%.1fs)', t(im1)));
end
plot(t(imQ), Vrel(imQ), 'd', 'Color', cMQ, 'MarkerSize', 9, ...
    'MarkerFaceColor', cMQ, 'DisplayName', sprintf('Max-Q (t=%.1fs)', t(imQ)));
xlabel('Time [s]'); ylabel('Velocity [m/s]');
title('Velocity Magnitude vs Time  (ND simulation)');
legend('Location', 'northwest'); grid on; hold off;

% ---- Figure 4: Mass vs Time ----
figure('Name', 'Mass [ND sim]', 'Position', [200 250 700 450]);
plot(t, mass/1e3, 'b', 'LineWidth', 1.5); hold on;
xline(t_end1, '--', 'Color', cP1, 'LineWidth', 1, 'HandleVisibility', 'off');
xline(t_end2, '--', 'Color', cP2, 'LineWidth', 1, 'HandleVisibility', 'off');
plot(t(ip1), mass(ip1)/1e3, 'v', 'Color', cP1, 'MarkerSize', 9, ...
    'MarkerFaceColor', cP1, 'DisplayName', 'End Vert. Ascent');
plot(t(ip2), mass(ip2)/1e3, 'v', 'Color', cP2, 'MarkerSize', 9, ...
    'MarkerFaceColor', cP2, 'DisplayName', 'End Pitchover');
xlabel('Time [s]'); ylabel('Mass [t]');
title('Vehicle Mass vs Time  (ND simulation)');
legend('Mass', 'Location', 'northeast'); grid on; hold off;

% ---- Figure 5: Dynamic Pressure vs Time ----
figure('Name', 'Dynamic Pressure [ND sim]', 'Position', [250 200 700 450]);
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
title('Dynamic Pressure vs Time  (ND simulation)');
legend('Location', 'northeast'); grid on; hold off;

% ---- Figure 6: Mach Number vs Time ----
figure('Name', 'Mach Number [ND sim]', 'Position', [300 150 700 450]);
plot(t, Mach, 'b', 'LineWidth', 1.5, 'DisplayName', 'Mach'); hold on;
yline(1, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
xline(t_end1, '--', 'Color', cP1, 'LineWidth', 1, 'HandleVisibility', 'off');
xline(t_end2, '--', 'Color', cP2, 'LineWidth', 1, 'HandleVisibility', 'off');
if ~isempty(im1)
    plot(t(im1), 1, '^', 'Color', cM1, 'MarkerSize', 9, ...
        'MarkerFaceColor', cM1, 'DisplayName', sprintf('Mach 1 (t=%.1fs)', t(im1)));
end
xlabel('Time [s]'); ylabel('Mach Number [-]');
title('Mach Number vs Time  (ND simulation)');
legend('Location', 'northwest'); grid on; hold off;

% ---- Figure 7: Angles vs Time ----
figure('Name', 'Angles [ND sim]', 'Position', [350 100 700 450]);
plot(t, gammaT, 'b',  'LineWidth', 1.5, 'DisplayName', '\gamma_T (Thrust)'); hold on;
plot(t, gammaA, 'r--','LineWidth', 1.5, 'DisplayName', '\gamma_A (Aerodynamic)');
xline(t_end1, '--', 'Color', cP1, 'LineWidth', 1, 'HandleVisibility', 'off');
xline(t_end2, '--', 'Color', cP2, 'LineWidth', 1, 'HandleVisibility', 'off');
xlabel('Time [s]'); ylabel('Angle [deg]');
title('Thrust and Aerodynamic Flight-Path Angles  (ND simulation)');
legend('Location', 'best'); grid on; hold off;

% ---- Figure 8: Ground Track ----
figure('Name', 'Ground Track [ND sim]', 'Position', [400 50 700 450]);
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
title('Ground Track  (ND simulation)');
legend('Trajectory', 'Location', 'best'); grid on; hold off;

% ---- Figure 9: Non-dimensional state variables ----
figure('Name', 'ND State', 'Position', [450 50 900 600]);

subplot(3,2,1);
plot(t_nd, Y_nd(:,1), 'b', 'LineWidth', 1.5);
xlabel('t^*'); ylabel('r^* = r/R_E');
title('Non-dim. radius'); grid on;

subplot(3,2,2);
plot(t_nd, Y_nd(:,4), 'b', 'LineWidth', 1.5); hold on;
plot(t_nd, Y_nd(:,5), 'r--', 'LineWidth', 1.5);
plot(t_nd, Y_nd(:,6), 'g:', 'LineWidth', 1.5);
xlabel('t^*'); ylabel('V^* = V/V_{ref}');
legend('u^*','v^*','w^*'); title('Non-dim. velocity components'); grid on;

subplot(3,2,3);
plot(t_nd, Y_nd(:,7), 'b', 'LineWidth', 1.5);
xlabel('t^*'); ylabel('m^* = m/m_0');
title('Non-dim. mass'); grid on;

subplot(3,2,4);
plot(t_nd, Y_nd(:,1) - 1, 'b', 'LineWidth', 1.5);
xlabel('t^*'); ylabel('h^* = h/R_E');
title('Non-dim. altitude'); grid on;

subplot(3,2,5);
Vmag_nd = sqrt(Y_nd(:,4).^2 + Y_nd(:,5).^2 + Y_nd(:,6).^2);
plot(t_nd, Vmag_nd, 'b', 'LineWidth', 1.5);
xlabel('t^*'); ylabel('|V^*|');
title('Non-dim. inertial speed'); grid on;

subplot(3,2,6);
% Visualise arc boundaries in nd time
yyaxis left;
plot(t_nd, Y_nd(:,7), 'b', 'LineWidth', 1.5);
ylabel('m^* [-]');
yyaxis right;
plot(t_nd, Vmag_nd, 'r', 'LineWidth', 1.5);
ylabel('|V^*| [-]');
xline(t1_nd, 'k--', 'LineWidth', 1);
xline(t2_nd, 'k--', 'LineWidth', 1);
xlabel('t^* (global nd time)');
title('Mass & speed with arc boundaries'); grid on;

sgtitle('Non-dimensional state variables  (L_{ref}=R_E,  V_{ref}=7800 m/s)');

%% ========================================================================
%  NON-DIMENSIONAL EQUATIONS OF MOTION  (local function)
%  ========================================================================

function dydt = eom_nd(tau, y, par)
% Non-dimensional EOM for a single flight arc.
%
% Inputs
%   tau  : arc-local dimensionless time in [0, 1]
%            tau = (t_dim - t_arc_start_dim) / arc_duration_dim
%   y    : non-dimensional state [r*, theta, phi, u*, v*, w*, m*]
%   par  : parameter structure (all quantities non-dimensional unless noted)
%
% The global non-dimensional time is
%   t_nd = par.t_arc_start + tau * par.Delta
% The corresponding dimensional time is
%   t_dim = t_nd * par.T_ref
%
% EOM scaled by arc duration:
%   dy*/dtau = par.Delta * (dy*/dt*)
% so that each arc spans tau in [0,1].

    % Global nd time (used for pitchover angle)
    t_nd = par.t_arc_start + tau * par.Delta;

    % Unpack non-dimensional state
    r   = y(1);   % r / L_ref
    phi = y(3);   % latitude [rad]  (dimensionless angle)
    u   = y(4);   % u / V_ref
    v   = y(5);   % v / V_ref
    w   = y(6);   % w / V_ref
    m   = y(7);   % m / m_ref

    % Non-dimensional altitude  (R_E* = 1 since L_ref = R_E)
    alt = r - 1;

    % Exponential atmosphere (non-dimensional)
    rho  = par.rho0 * exp(-alt / par.H);   % rho / rho_ref
    patm = par.p0   * exp(-alt / par.H);   % p / p_ref

    % Atmospheric co-rotation velocity (non-dimensional)
    v_atm = par.omegaE * r * cos(phi);     % v_atm / V_ref

    % Relative velocity components
    ur = u;
    vr = v - v_atm;
    wr = w;
    Vrel = sqrt(ur^2 + vr^2 + wr^2);      % |V_rel| / V_ref

    % Gravity  (non-dimensional, radial only)
    %   g_nd = -(mu / r^2) / (V_ref^2/L_ref)  = -mu_nd / r_nd^2
    g_u = -par.mu / r^2;

    % Aerodynamic drag (non-dimensional)
    %   D_nd/m_nd = -0.5 * rho_nd * Sref_nd * CD * Vrel_nd * [ur,vr,wr]_nd / m_nd
    if Vrel > 1e-10
        Dc = -0.5 * rho * par.Sref * par.CD * Vrel;
        Du = Dc * ur;
        Dv = Dc * vr;
        Dw = Dc * wr;
    else
        Du = 0;  Dv = 0;  Dw = 0;
    end

    % Thrust magnitude (non-dimensional)
    %   Tmag_nd = Tvac_nd - patm_nd * Aex_nd
    Tmag = par.Tvac - patm * par.Aex;

    % Thrust direction (depends on current arc/phase)
    switch par.phase

        case 1
            % Vertical ascent: thrust straight up (radial)
            Tu = Tmag;  Tv = 0;  Tw = 0;

        case 2
            % Pitchover: linear elevation 90 -> 89.5 deg over 10 s
            % Elevation angle defined in dimensional time
            t_dim = t_nd * par.T_ref;
            gT    = deg2rad(90 - 0.05 * (t_dim - par.t_end1));
            psi   = deg2rad(90);          % azimuth = East
            Tu    = Tmag * sin(gT);
            Tv    = Tmag * cos(gT) * sin(psi);
            Tw    = Tmag * cos(gT) * cos(psi);

        case 3
            % Zero-lift gravity turn: thrust along relative velocity
            if Vrel > 1e-10
                Tu = Tmag * ur / Vrel;
                Tv = Tmag * vr / Vrel;
                Tw = Tmag * wr / Vrel;
            else
                Tu = Tmag;  Tv = 0;  Tw = 0;
            end
    end

    % --- Kinematic equations (non-dimensional, same structure as dim.) ---
    dr     = u;
    dtheta = v / (r * cos(phi));
    dphi   = w / r;

    % --- Dynamic equations (transport terms + non-dimensional forces/m) ---
    du = (v^2 + w^2) / r  +  g_u  +  (Tu + Du) / m;
    dv = (-u*v + v*w*tan(phi)) / r  +  (Tv + Dv) / m;
    dw = (-u*w - v^2*tan(phi)) / r  +  (Tw + Dw) / m;

    % --- Mass equation ---
    dm = -par.Qdot;

    % Scale by arc duration: d()/dtau = Delta * d()/dt_nd
    dydt = par.Delta * [dr; dtheta; dphi; du; dv; dw; dm];
end
