%% main2.m - Falcon 9 First Stage Trajectory Simulation (Non-Dimensional)
%  Classwork n. 1 - Simulation of Falcon 9 First Stage Motion
%  3-DoF point-mass model in spherical coordinates
%  Velocity in Up-East-North (UEN) frame
%
%  Non-dimensionalization:
%    L_ref = RE = 6378137 m        (Earth equatorial radius)
%    V_ref      = 7800 m/s         (first cosmic velocity, ~7.8 km/s)
%    T_ref = L_ref / V_ref         (time reference)
%    m_ref = m0                    (initial total mass)
%
%  Derived reference scales (implicit):
%    a_ref = V_ref^2 / L_ref       (acceleration)
%    F_ref = m_ref * a_ref         (force)
%    p_ref = F_ref / L_ref^2       (pressure)
%
%  Integration variable:
%    A single parameter tau in [0, 3] covers all three arcs.
%    Each arc occupies one unit of tau so that within each arc
%    the arc-local variable is in [0, 1]:
%
%      tau in [0,1]  -> Arc 1 (vertical ascent):
%                       t* =  tau          * Delta1
%      tau in [1,2]  -> Arc 2 (pitchover):
%                       t* = (tau - 1)     * Delta2  +  t1*
%      tau in [2,3]  -> Arc 3 (gravity turn):
%                       t* = (tau - 2)     * Delta3  +  t2*
%
%    where t* = t_dim / T_ref  (global non-dimensional time)
%    and   Delta_k = arc_duration_dim / T_ref.
%
%    The EOM is scaled accordingly:  dy*/dtau = Delta_k * (dy*/dt*)
%    ode45 integrates once over tspan = [0, 3].

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
a_sound   = sqrt(gamma_air * Rgas * Tamb);   % Speed of sound [m/s]

% --- First Stage ---
mdry1  = 22200;              % Dry mass [kg]
mp1    = 410900;             % Propellant mass [kg]
Tvac1  = 8227e3;             % Vacuum thrust [N]
tb1    = 162;                % Burn time [s]
cvac1  = 3244;               % Vacuum exhaust velocity [m/s]
Aex1   = 11.039;             % Nozzle exit area [m^2]
Qdot1  = Tvac1 / cvac1;     % Mass flow rate [kg/s]

% --- Second Stage (mass contribution only) ---
mdry2 = 4000;
mp2   = 107500;

% --- Aerodynamics ---
CD   = 0.329;
Sref = 10.52;                % Reference area [m^2]

% --- Payload and fairing ---
mfair = 1700;
mpay  = 22800;

% --- Phase timing ---
t_end1 = 5;                  % End of vertical ascent [s]
t_end2 = 15;                 % End of pitchover [s]

%% ========================================================================
%  LAUNCH SITE — Kennedy Space Center
%  ========================================================================

lat0 = deg2rad(28.573469);
lon0 = deg2rad(-80.651070);

%% ========================================================================
%  REFERENCE QUANTITIES
%  ========================================================================

L_ref = RE;            % [m]
V_ref = 7800;          % [m/s]
T_ref = L_ref / V_ref; % [s]

m0    = mdry1 + mp1 + mdry2 + mp2 + mfair + mpay;
m_ref = m0;            % [kg]

fprintf('=== Reference Scales ===\n');
fprintf('  L_ref = %.4f km\n', L_ref/1e3);
fprintf('  V_ref = %.1f  m/s\n', V_ref);
fprintf('  T_ref = %.4f s\n',   T_ref);
fprintf('  m_ref = %.0f  kg\n', m_ref);
fprintf('========================\n\n');

%% ========================================================================
%  NON-DIMENSIONAL PARAMETERS
%  ========================================================================

mu_nd     = mu     / (V_ref^2 * L_ref);          % gravitational parameter
omegaE_nd = omegaE * T_ref;                       % Earth angular velocity
rho0_nd   = rho0   * L_ref^3 / m_ref;            % sea-level density
p0_nd     = p0     * L_ref^3 / (m_ref * V_ref^2); % sea-level pressure
H_nd      = Hscale / L_ref;                       % scale height
Tvac_nd   = Tvac1  * L_ref  / (m_ref * V_ref^2); % vacuum thrust
Aex_nd    = Aex1   / L_ref^2;                    % nozzle exit area
Qdot_nd   = Qdot1  * T_ref  / m_ref;             % mass flow rate
Sref_nd   = Sref   / L_ref^2;                    % reference area
% CD is already dimensionless

%% ========================================================================
%  NON-DIMENSIONAL INITIAL CONDITIONS
%  ========================================================================

r0_nd = 1;                                         % RE / L_ref = 1
u0_nd = 0;
v0_nd = omegaE * RE * cos(lat0) / V_ref;           % Earth-rotation velocity
w0_nd = 0;
m0_nd = 1;                                         % m0 / m_ref = 1

y0_nd = [r0_nd; lon0; lat0; u0_nd; v0_nd; w0_nd; m0_nd];

%% ========================================================================
%  ARC DURATIONS IN NON-DIMENSIONAL TIME (Delta_k = duration / T_ref)
%  ========================================================================

t1_nd  = t_end1 / T_ref;    % global nd time at end of arc 1
t2_nd  = t_end2 / T_ref;    % global nd time at end of arc 2
tb_nd  = tb1    / T_ref;    % global nd time at MECO

Delta1 = t1_nd;              % arc 1 nd duration
Delta2 = t2_nd - t1_nd;     % arc 2 nd duration
Delta3 = tb_nd - t2_nd;     % arc 3 nd duration

%% ========================================================================
%  PARAMETER STRUCTURE
%  ========================================================================

par.mu      = mu_nd;
par.omegaE  = omegaE_nd;
par.rho0    = rho0_nd;
par.p0      = p0_nd;
par.H       = H_nd;
par.Tvac    = Tvac_nd;
par.Aex     = Aex_nd;
par.Qdot    = Qdot_nd;
par.CD      = CD;
par.Sref    = Sref_nd;
par.Delta1  = Delta1;
par.Delta2  = Delta2;
par.Delta3  = Delta3;
par.t1_nd   = t1_nd;        % global nd boundary between arcs 1 and 2
par.t2_nd   = t2_nd;        % global nd boundary between arcs 2 and 3
par.T_ref   = T_ref;        % [s]  to recover dim. time for pitchover law
par.t_end1  = t_end1;       % [s]  pitchover start

%% ========================================================================
%  SINGLE ode45 CALL — tau in [0, 3]
%
%  tau in [0,1]: Arc 1  ->  dy*/dtau = Delta1 * f_nd(t*, y*)
%  tau in [1,2]: Arc 2  ->  dy*/dtau = Delta2 * f_nd(t*, y*)
%  tau in [2,3]: Arc 3  ->  dy*/dtau = Delta3 * f_nd(t*, y*)
%
%  Arc boundaries are exact (Events function stops the step there).
%  ========================================================================

opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12, ...
              'Events', @arc_boundary_events);

[tau_sol, Y_nd] = ode45(@(tau, y) eom_nd(tau, y, par), [0 3], y0_nd, opts);

%% ========================================================================
%  RECOVER GLOBAL ND TIME AND DIMENSIONAL QUANTITIES
%  ========================================================================

% Global nd time from tau
t_nd = zeros(size(tau_sol));
mask1 =              tau_sol <= 1;
mask2 = tau_sol > 1 & tau_sol <= 2;
mask3 = tau_sol > 2;
t_nd(mask1) =  tau_sol(mask1)       * Delta1;
t_nd(mask2) = (tau_sol(mask2) - 1)  * Delta2  + t1_nd;
t_nd(mask3) = (tau_sol(mask3) - 2)  * Delta3  + t2_nd;

% Dimensional time [s]
t = t_nd * T_ref;

% State variables (dimensional)
r     = Y_nd(:,1) * L_ref;
theta = Y_nd(:,2);
phi   = Y_nd(:,3);
uvel  = Y_nd(:,4) * V_ref;
vvel  = Y_nd(:,5) * V_ref;
wvel  = Y_nd(:,6) * V_ref;
mass  = Y_nd(:,7) * m_ref;

% Altitude
h = r - RE;

% Inertial and relative velocity
Vmag = sqrt(uvel.^2 + vvel.^2 + wvel.^2);
urel = uvel;
vrel = vvel - omegaE .* r .* cos(phi);
wrel = wvel;
Vrel = sqrt(urel.^2 + vrel.^2 + wrel.^2);

% Atmospheric quantities along trajectory
rho_traj = rho0 * exp(-h / Hscale);
p_traj   = p0   * exp(-h / Hscale);
qdyn     = 0.5 * rho_traj .* Vrel.^2;
Mach     = Vrel / a_sound;

% Thrust and aerodynamic angles
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

% Mission events
[~, ip1] = min(abs(t - t_end1));
[~, ip2] = min(abs(t - t_end2));
im1      = find(Mach >= 1, 1, 'first');
[qmax, imQ] = max(qdyn);

% ECI -> local ENU
Xeci = r .* cos(phi) .* cos(theta);
Yeci = r .* cos(phi) .* sin(theta);
Zeci = r .* sin(phi);
X0   = RE * cos(lat0) * cos(lon0);
Y0   = RE * cos(lat0) * sin(lon0);
Z0   = RE * sin(lat0);
R_enu = [-sin(lon0),            cos(lon0),             0;
         -sin(lat0)*cos(lon0), -sin(lat0)*sin(lon0),   cos(lat0);
          cos(lat0)*cos(lon0),  cos(lat0)*sin(lon0),   sin(lat0)];
dR_enu   = (R_enu * [Xeci-X0, Yeci-Y0, Zeci-Z0]')';
East_km  = dR_enu(:,1)/1e3;
North_km = dR_enu(:,2)/1e3;
Up_km    = dR_enu(:,3)/1e3;

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
fprintf('    r* = %.8f    (r / R_E)\n',   Y_nd(end,1));
fprintf('    u* = %.8f    (u / V_ref)\n', Y_nd(end,4));
fprintf('    v* = %.8f    (v / V_ref)\n', Y_nd(end,5));
fprintf('    w* = %.8f    (w / V_ref)\n', Y_nd(end,6));
fprintf('    m* = %.8f    (m / m0)\n',    Y_nd(end,7));
fprintf('    t* = %.8f    (t / T_ref)\n', t_nd(end));
fprintf('====================================================\n');

%% ========================================================================
%  PLOTS
%  ========================================================================

cP1 = [0.2 0.7 0.2];
cP2 = [0.9 0.6 0.0];
cM1 = [0.85 0.0 0.0];
cMQ = [0.0 0.2 0.8];

% ---- Figure 1: 3D Trajectory ----
figure('Name','3D Trajectory [ND]','Position',[50 400 700 550]);
plot3(East_km,North_km,Up_km,'b','LineWidth',1.8); hold on;
plot3(East_km(1),  North_km(1),  Up_km(1),  'go','MarkerSize',12,'MarkerFaceColor','g');
plot3(East_km(end),North_km(end),Up_km(end),'rs','MarkerSize',12,'MarkerFaceColor','r');
plot3(East_km(ip1),North_km(ip1),Up_km(ip1),'v','Color',cP1,'MarkerSize',10,'MarkerFaceColor',cP1);
plot3(East_km(ip2),North_km(ip2),Up_km(ip2),'v','Color',cP2,'MarkerSize',10,'MarkerFaceColor',cP2);
if ~isempty(im1)
    plot3(East_km(im1),North_km(im1),Up_km(im1),'^','Color',cM1,'MarkerSize',10,'MarkerFaceColor',cM1);
end
plot3(East_km(imQ),North_km(imQ),Up_km(imQ),'d','Color',cMQ,'MarkerSize',10,'MarkerFaceColor',cMQ);
xlabel('East [km]'); ylabel('North [km]'); zlabel('Up [km]');
title('Falcon 9 First Stage – 3D Trajectory  (ND simulation)');
legend('Trajectory','Launch','MECO','End Vert. Ascent','End Pitchover', ...
       'Mach 1','Max-Q','Location','best');
grid on; view([-35 25]); hold off;

% ---- Figure 2: Altitude vs Time ----
figure('Name','Altitude [ND]','Position',[100 350 700 450]);
plot(t,h/1e3,'b','LineWidth',1.5); hold on;
xline(t_end1,'--','Color',cP1,'LineWidth',1,'HandleVisibility','off');
xline(t_end2,'--','Color',cP2,'LineWidth',1,'HandleVisibility','off');
plot(t(ip1),h(ip1)/1e3,'v','Color',cP1,'MarkerSize',9,'MarkerFaceColor',cP1,'DisplayName','End Vert. Ascent');
plot(t(ip2),h(ip2)/1e3,'v','Color',cP2,'MarkerSize',9,'MarkerFaceColor',cP2,'DisplayName','End Pitchover');
if ~isempty(im1)
    plot(t(im1),h(im1)/1e3,'^','Color',cM1,'MarkerSize',9,'MarkerFaceColor',cM1,'DisplayName',sprintf('Mach 1 (t=%.1fs)',t(im1)));
end
plot(t(imQ),h(imQ)/1e3,'d','Color',cMQ,'MarkerSize',9,'MarkerFaceColor',cMQ,'DisplayName',sprintf('Max-Q (t=%.1fs)',t(imQ)));
xlabel('Time [s]'); ylabel('Altitude [km]');
title('Altitude vs Time  (ND simulation)');
legend('Altitude','Location','northwest'); grid on; hold off;

% ---- Figure 3: Velocity ----
figure('Name','Velocity [ND]','Position',[150 300 700 450]);
plot(t,Vmag,'b','LineWidth',1.5,'DisplayName','Inertial |V|'); hold on;
plot(t,Vrel,'r--','LineWidth',1.5,'DisplayName','Relative |V_{rel}|');
xline(t_end1,'--','Color',cP1,'LineWidth',1,'HandleVisibility','off');
xline(t_end2,'--','Color',cP2,'LineWidth',1,'HandleVisibility','off');
if ~isempty(im1)
    plot(t(im1),Vrel(im1),'^','Color',cM1,'MarkerSize',9,'MarkerFaceColor',cM1,'DisplayName',sprintf('Mach 1 (t=%.1fs)',t(im1)));
end
plot(t(imQ),Vrel(imQ),'d','Color',cMQ,'MarkerSize',9,'MarkerFaceColor',cMQ,'DisplayName',sprintf('Max-Q (t=%.1fs)',t(imQ)));
xlabel('Time [s]'); ylabel('Velocity [m/s]');
title('Velocity Magnitude vs Time  (ND simulation)');
legend('Location','northwest'); grid on; hold off;

% ---- Figure 4: Mass ----
figure('Name','Mass [ND]','Position',[200 250 700 450]);
plot(t,mass/1e3,'b','LineWidth',1.5); hold on;
xline(t_end1,'--','Color',cP1,'LineWidth',1,'HandleVisibility','off');
xline(t_end2,'--','Color',cP2,'LineWidth',1,'HandleVisibility','off');
plot(t(ip1),mass(ip1)/1e3,'v','Color',cP1,'MarkerSize',9,'MarkerFaceColor',cP1,'DisplayName','End Vert. Ascent');
plot(t(ip2),mass(ip2)/1e3,'v','Color',cP2,'MarkerSize',9,'MarkerFaceColor',cP2,'DisplayName','End Pitchover');
xlabel('Time [s]'); ylabel('Mass [t]');
title('Vehicle Mass vs Time  (ND simulation)');
legend('Mass','Location','northeast'); grid on; hold off;

% ---- Figure 5: Dynamic Pressure ----
figure('Name','Dynamic Pressure [ND]','Position',[250 200 700 450]);
plot(t,qdyn/1e3,'b','LineWidth',1.5,'DisplayName','q_{dyn}'); hold on;
xline(t_end1,'--','Color',cP1,'LineWidth',1,'HandleVisibility','off');
xline(t_end2,'--','Color',cP2,'LineWidth',1,'HandleVisibility','off');
if ~isempty(im1)
    plot(t(im1),qdyn(im1)/1e3,'^','Color',cM1,'MarkerSize',9,'MarkerFaceColor',cM1,'DisplayName',sprintf('Mach 1 (t=%.1fs)',t(im1)));
end
plot(t(imQ),qmax/1e3,'d','Color',cMQ,'MarkerSize',12,'MarkerFaceColor',cMQ,'DisplayName',sprintf('Max-Q = %.1f kPa',qmax/1e3));
xlabel('Time [s]'); ylabel('Dynamic Pressure [kPa]');
title('Dynamic Pressure vs Time  (ND simulation)');
legend('Location','northeast'); grid on; hold off;

% ---- Figure 6: Mach Number ----
figure('Name','Mach [ND]','Position',[300 150 700 450]);
plot(t,Mach,'b','LineWidth',1.5,'DisplayName','Mach'); hold on;
yline(1,'k--','LineWidth',1,'HandleVisibility','off');
xline(t_end1,'--','Color',cP1,'LineWidth',1,'HandleVisibility','off');
xline(t_end2,'--','Color',cP2,'LineWidth',1,'HandleVisibility','off');
if ~isempty(im1)
    plot(t(im1),1,'^','Color',cM1,'MarkerSize',9,'MarkerFaceColor',cM1,'DisplayName',sprintf('Mach 1 (t=%.1fs)',t(im1)));
end
xlabel('Time [s]'); ylabel('Mach Number [-]');
title('Mach Number vs Time  (ND simulation)');
legend('Location','northwest'); grid on; hold off;

% ---- Figure 7: Angles ----
figure('Name','Angles [ND]','Position',[350 100 700 450]);
plot(t,gammaT,'b','LineWidth',1.5,'DisplayName','\gamma_T (Thrust)'); hold on;
plot(t,gammaA,'r--','LineWidth',1.5,'DisplayName','\gamma_A (Aerodynamic)');
xline(t_end1,'--','Color',cP1,'LineWidth',1,'HandleVisibility','off');
xline(t_end2,'--','Color',cP2,'LineWidth',1,'HandleVisibility','off');
xlabel('Time [s]'); ylabel('Angle [deg]');
title('Thrust and Flight-Path Angles  (ND simulation)');
legend('Location','best'); grid on; hold off;

% ---- Figure 8: Ground Track ----
figure('Name','Ground Track [ND]','Position',[400 50 700 450]);
plot(lon_ground,lat_ground,'b','LineWidth',2); hold on;
plot(lon_ground(1),lat_ground(1),'go','MarkerSize',12,'MarkerFaceColor','g','DisplayName','Launch (KSC)');
plot(lon_ground(end),lat_ground(end),'rs','MarkerSize',12,'MarkerFaceColor','r','DisplayName','MECO');
if ~isempty(im1)
    plot(lon_ground(im1),lat_ground(im1),'^','Color',cM1,'MarkerSize',9,'MarkerFaceColor',cM1,'DisplayName','Mach 1');
end
plot(lon_ground(imQ),lat_ground(imQ),'d','Color',cMQ,'MarkerSize',9,'MarkerFaceColor',cMQ,'DisplayName','Max-Q');
xlabel('Longitude [deg]'); ylabel('Latitude [deg]');
title('Ground Track  (ND simulation)');
legend('Trajectory','Location','best'); grid on; hold off;

% ---- Figure 9: Non-dimensional state variables ----
figure('Name','ND State','Position',[450 30 900 580]);
tlo = tiledlayout(3, 2);

nexttile;
plot(tau_sol, Y_nd(:,1), 'b', 'LineWidth',1.5);
xline(1,'k--'); xline(2,'k--');
xlabel('\tau'); ylabel('r^* = r / R_E');
title('Non-dim. radius'); grid on;

nexttile;
plot(tau_sol, Y_nd(:,4),'b','LineWidth',1.5); hold on;
plot(tau_sol, Y_nd(:,5),'r--','LineWidth',1.5);
plot(tau_sol, Y_nd(:,6),'g:','LineWidth',1.5);
xline(1,'k--'); xline(2,'k--');
xlabel('\tau'); ylabel('V^* = V / V_{ref}');
legend('u^*','v^*','w^*'); title('Non-dim. velocity components'); grid on;

nexttile;
plot(tau_sol, Y_nd(:,7), 'b', 'LineWidth',1.5);
xline(1,'k--'); xline(2,'k--');
xlabel('\tau'); ylabel('m^* = m / m_0');
title('Non-dim. mass'); grid on;

nexttile;
plot(tau_sol, Y_nd(:,1)-1, 'b', 'LineWidth',1.5);
xline(1,'k--'); xline(2,'k--');
xlabel('\tau'); ylabel('h^* = h / R_E');
title('Non-dim. altitude'); grid on;

nexttile;
Vmag_nd = sqrt(Y_nd(:,4).^2 + Y_nd(:,5).^2 + Y_nd(:,6).^2);
plot(tau_sol, Vmag_nd, 'b', 'LineWidth',1.5);
xline(1,'k--','Arc 1|2'); xline(2,'k--','Arc 2|3');
xlabel('\tau  (0\rightarrow1: arc 1,  1\rightarrow2: arc 2,  2\rightarrow3: arc 3)');
ylabel('|V^*|'); title('Non-dim. inertial speed'); grid on;

nexttile;
plot(tau_sol, t, 'b', 'LineWidth',1.5);
xline(1,'k--'); xline(2,'k--');
xlabel('\tau'); ylabel('t_{dim}  [s]');
title('Dimensional time vs \tau  (shows per-arc scaling)'); grid on;

title(tlo, sprintf( ...
    'Non-dimensional state   (L_{ref}=R_E,  V_{ref}=%g m/s,  T_{ref}=%.1f s)', ...
    V_ref, T_ref));

%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

function dydt = eom_nd(tau, y, par)
% Non-dimensional EOM for a single continuous ode45 call over tau in [0,3].
%
% tau in [0,1] -> Arc 1 (vertical ascent),  dy*/dtau = Delta1 * f_nd
% tau in [1,2] -> Arc 2 (pitchover),        dy*/dtau = Delta2 * f_nd
% tau in [2,3] -> Arc 3 (gravity turn),     dy*/dtau = Delta3 * f_nd
%
% Within each arc the local variable is tau_k = tau - (k-1) in [0,1].
% The global nd time is recovered as:
%   t* = t_arc_start* + tau_k * Delta_k

    % --- Identify current arc and recover t* ---
    if tau <= 1
        Delta       = par.Delta1;
        t_arc_start = 0;
        tau_k       = tau;
        phase       = 1;
    elseif tau <= 2
        Delta       = par.Delta2;
        t_arc_start = par.t1_nd;
        tau_k       = tau - 1;
        phase       = 2;
    else
        Delta       = par.Delta3;
        t_arc_start = par.t2_nd;
        tau_k       = tau - 2;
        phase       = 3;
    end
    t_nd = t_arc_start + tau_k * Delta;

    % --- Unpack non-dimensional state ---
    r   = y(1);   % r / L_ref
    phi = y(3);   % latitude [rad]
    u   = y(4);   % u / V_ref
    v   = y(5);   % v / V_ref
    w   = y(6);   % w / V_ref
    m   = y(7);   % m / m_ref

    % Non-dimensional altitude  (R_E* = 1 since L_ref = R_E)
    alt = r - 1;

    % Atmosphere
    rho  = par.rho0 * exp(-alt / par.H);
    patm = par.p0   * exp(-alt / par.H);

    % Atmospheric co-rotation velocity (nd)
    v_atm = par.omegaE * r * cos(phi);

    % Relative velocity components (nd)
    ur   = u;
    vr   = v - v_atm;
    wr   = w;
    Vrel = sqrt(ur^2 + vr^2 + wr^2);

    % Gravity (nd, radial)
    g_u = -par.mu / r^2;

    % Aerodynamic drag (nd)
    %   D_nd/m_nd = -0.5 * rho_nd * Sref_nd * CD * Vrel_nd * [ur,vr,wr]_nd / m_nd
    if Vrel > 1e-10
        Dc = -0.5 * rho * par.Sref * par.CD * Vrel;
        Du = Dc * ur;   Dv = Dc * vr;   Dw = Dc * wr;
    else
        Du = 0;   Dv = 0;   Dw = 0;
    end

    % Thrust magnitude (nd):  Tmag* = Tvac* - p* * Aex*
    Tmag = par.Tvac - patm * par.Aex;

    % Thrust direction
    switch phase
        case 1
            % Vertical ascent — thrust radially upward
            Tu = Tmag;   Tv = 0;   Tw = 0;

        case 2
            % Pitchover — linear elevation decrease in dimensional time
            t_dim = t_nd * par.T_ref;
            gT    = deg2rad(90 - 0.05 * (t_dim - par.t_end1));
            % Azimuth psi = 90 deg (East): sin(psi)=1, cos(psi)=0
            Tu    = Tmag * sin(gT);
            Tv    = Tmag * cos(gT);
            Tw    = 0;

        case 3
            % Gravity turn — thrust along relative velocity
            if Vrel > 1e-10
                Tu = Tmag * ur / Vrel;
                Tv = Tmag * vr / Vrel;
                Tw = Tmag * wr / Vrel;
            else
                Tu = Tmag;   Tv = 0;   Tw = 0;
            end
    end

    % Kinematic equations (nd, same structure as dimensional)
    dr     = u;
    dtheta = v / (r * cos(phi));
    dphi   = w / r;

    % Dynamic equations (transport terms + nd force / nd mass)
    du = (v^2 + w^2) / r  +  g_u  +  (Tu + Du) / m;
    dv = (-u*v + v*w*tan(phi)) / r  +  (Tv + Dv) / m;
    dw = (-u*w - v^2*tan(phi)) / r  +  (Tw + Dw) / m;

    % Mass equation
    dm = -par.Qdot;

    % Chain rule: dy*/dtau = Delta_k * (dy*/dt*)
    dydt = Delta * [dr; dtheta; dphi; du; dv; dw; dm];
end

% --------------------------------------------------------------------------

function [value, isterminal, direction] = arc_boundary_events(tau, ~, ~)
% Locate tau = 1 and tau = 2 exactly so ode45 never straddles an arc
% boundary (where the scaling Delta changes discontinuously).
% isterminal = 0  -> solver continues past the event.
    value      = [tau - 1;   tau - 2];
    isterminal = [0;          0];
    direction  = [+1;         +1];
end
