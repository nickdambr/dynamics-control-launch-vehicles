%% TASK 3 - Vertical climb + burn + coasting arc
%  Mission sequence: vertical - burn - coast
%  Search for optimal engine cutoff time to maximize payload

clear; close all; clc;

%% Parameters
c   = 0.6;
eta = 0.1;
y1  = 0.0001;   % vertical climb altitude
yf  = 0.04;     % target orbit altitude
Q   = 2;         % mass flow rate
T   = c * Q;

opts_ode = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
opts_fs  = optimoptions('fsolve', 'Display', 'iter', ...
    'MaxIterations', 500, 'MaxFunctionEvaluations', 5000, ...
    'FunctionTolerance', 1e-10, 'StepTolerance', 1e-10);

%% Phase 1: Vertical climb to y1
fprintf('Phase 1: Vertical climb to y1 = %.4f\n', y1);
ic_vert = [0; 0; 1];
opts_vert = odeset('RelTol', 1e-12, 'AbsTol', 1e-14, ...
    'Events', @(t,z) event_altitude(t, z, y1));
[T_vert, Z_vert] = ode45(@(t,z) ode_vertical(t, z, T, Q), [0 1], ic_vert, opts_vert);

t1  = T_vert(end);
vy1 = Z_vert(end,2);
m1  = Z_vert(end,3);

fprintf('  t1 = %.6f,  vy1 = %.6f,  m1 = %.6f\n', t1, vy1, m1);

%% Phase 2+3: Burn + Coast
% BVP unknowns: [lam_vx0, lam_vy0, lam_y, lam_m0, t_burn]
%
% At engine cutoff (end of burn, time t_burn from start of Phase 2):
%   vx(tc) = 1               (coast keeps vx constant)
%   y(tc) + 0.5*vy(tc)^2 = yf  (coast ballistic condition)
%   lam_m(tc) = 1            (lam_m constant during coast, transversality)
%   lam_vy(tc) = lam_y*vy(tc)  (coast optimality from H_coast = 0)
%   S(tc) = |lam_v|/m - 1/c = 0  (switching function)

p.c  = c;
p.Q  = Q;
p.T  = T;
p.yf = yf;
p.x0  = 0;
p.y0  = y1;
p.vx0 = 0;
p.vy0 = vy1;
p.m0  = m1;

% Initial guess - start from Task 2 solution if available, otherwise manual
z_guess = [0.20; 1.2; 4.0; -0.5; 0.28];

fprintf('\nSolving BVP for burn + coast...\n');
[z_sol, fval, ef] = fsolve(@(z) shooting3(z, p, opts_ode), z_guess, opts_fs);

if ef <= 0
    fprintf('Trying alternative initial guesses...\n');
    guesses = {[0.15; 1.0; 3.0; -0.3; 0.25], ...
               [0.25; 1.5; 5.0; -0.7; 0.32], ...
               [0.18; 1.1; 3.5; -0.4; 0.27]};
    for gg = 1:length(guesses)
        [z_sol, fval, ef] = fsolve(@(z) shooting3(z, p, opts_ode), guesses{gg}, opts_fs);
        if ef > 0, break; end
    end
end

if ef > 0
    pp = p;
    pp.lam_vx0 = z_sol(1);
    pp.lam_vy0 = z_sol(2);
    pp.lam_y   = z_sol(3);
    t_burn = z_sol(5);

    % Re-integrate burn phase
    ic2 = [p.x0; p.y0; p.vx0; p.vy0; p.m0; z_sol(4)];
    [T2, Z2] = ode45(@(t,z) ode_burn(t, z, pp), linspace(0, t_burn, 500), ic2, opts_ode);

    % State at engine cutoff
    xc  = Z2(end,1);
    yc  = Z2(end,2);
    vxc = Z2(end,3);
    vyc = Z2(end,4);
    mc  = Z2(end,5);

    % Coast phase (analytical)
    t_coast = vyc;               % time to reach vy = 0
    tf_total = t1 + t_burn + t_coast;

    % Coast trajectory
    t_c = linspace(0, t_coast, 200)';
    x_coast = xc + vxc * t_c;
    y_coast = yc + vyc * t_c - 0.5 * t_c.^2;
    vx_coast = vxc * ones(size(t_c));
    vy_coast = vyc - t_c;

    mf_task3 = mc;  % mass doesn't change during coast

    fprintf('\n=== Task 3 Results (Q = %.1f, yf = %.2f) ===\n', Q, yf);
    fprintf('  Vertical climb time: %.6f\n', t1);
    fprintf('  Burn time:           %.6f\n', t_burn);
    fprintf('  Coast time:          %.6f\n', t_coast);
    fprintf('  Total flight time:   %.6f\n', tf_total);
    fprintf('  Final mass:          %.6f\n', mf_task3);
    fprintf('  Final altitude:      %.6f (target: %.4f)\n', y_coast(end), yf);
    fprintf('  Final vx:            %.6f (target: 1)\n', vx_coast(end));
    fprintf('  Final vy:            %.6f (target: 0)\n', vy_coast(end));
    fprintf('  Payload (approx):    %.6f\n', mf_task3*(1+eta) - eta);

    %% Plots
    figure('Name','Task 3 - Trajectory');
    hold on; grid on;
    % Vertical
    plot(zeros(size(T_vert)), Z_vert(:,1), 'r-', 'LineWidth', 2, 'DisplayName', 'Vertical');
    % Burn
    plot(Z2(:,1), Z2(:,2), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Burn');
    % Coast
    plot(x_coast, y_coast, 'g--', 'LineWidth', 1.5, 'DisplayName', 'Coast');
    xlabel('x (nondim)'); ylabel('y (nondim)');
    title('Task 3: Trajectory (vertical + burn + coast)');
    legend('Location','best');

    % Angles during burn
    phi_burn = zeros(length(T2),1);
    psi_burn = zeros(length(T2),1);
    for kk = 1:length(T2)
        lam_vy_k = z_sol(2) - z_sol(3) * T2(kk);
        phi_burn(kk) = atan2(lam_vy_k, z_sol(1));
        V_k = sqrt(Z2(kk,3)^2 + Z2(kk,4)^2);
        if V_k > 1e-10
            psi_burn(kk) = atan2(Z2(kk,4), Z2(kk,3));
        else
            psi_burn(kk) = phi_burn(kk);
        end
    end

    % Angles during coast
    psi_coast = atan2(vy_coast, vx_coast);

    figure('Name','Task 3 - Angles');
    hold on; grid on;
    % Vertical phase
    plot([0 t1], [90 90], 'b-', 'LineWidth', 1, 'HandleVisibility','off');
    % Burn phase
    plot(t1 + T2, rad2deg(phi_burn), 'b-', 'LineWidth', 1.5, 'DisplayName', '\phi (thrust)');
    plot(t1 + T2, rad2deg(psi_burn), 'r--', 'LineWidth', 1.5, 'DisplayName', '\psi (velocity)');
    % Coast phase (no thrust)
    plot(t1 + t_burn + t_c, rad2deg(psi_coast), 'r--', 'LineWidth', 1, 'HandleVisibility','off');
    xline(t1, '--k', 'End vertical', 'LabelOrientation', 'horizontal');
    xline(t1 + t_burn, '--k', 'Engine cutoff', 'LabelOrientation', 'horizontal');
    xlabel('Time (nondim)'); ylabel('Angle (deg)');
    title('Task 3: Angles');
    legend('Location','best');

    % Mass profile
    figure('Name','Task 3 - Mass');
    hold on; grid on;
    plot(T_vert, Z_vert(:,3), 'r-', 'LineWidth', 1.5, 'DisplayName','Vertical');
    plot(t1 + T2, Z2(:,5), 'b-', 'LineWidth', 1.5, 'DisplayName','Burn');
    plot(t1 + t_burn + t_c, mc*ones(size(t_c)), 'g--', 'LineWidth', 1.5, 'DisplayName','Coast');
    xlabel('Time (nondim)'); ylabel('Mass (nondim)');
    title('Task 3: Mass profile');
    legend('Location','best');

else
    fprintf('ERROR: BVP did not converge for Task 3.\n');
end

%% ===================== LOCAL FUNCTIONS =====================

function dz = ode_vertical(t, z, T, Q)
    vy = z(2); m = z(3);
    dz = [vy; T/m - 1; -Q];
end

function [value, isterminal, direction] = event_altitude(t, z, y_target)
    value = z(1) - y_target;
    isterminal = 1;
    direction = 1;
end

function res = shooting3(z0, p, opts_ode)
% Shooting residual for burn + coast
%   z0 = [lam_vx0; lam_vy0; lam_y; lam_m0; t_burn]
%
% Conditions at engine cutoff (tc):
%   1. vxc = 1
%   2. yc + 0.5*vyc^2 = yf
%   3. lam_m(tc) = 1
%   4. lam_vy(tc) = lam_y * vyc  (coast optimality)
%   5. S(tc) = |lam_v|/mc - 1/c = 0  (switching function)

    lam_vx0 = z0(1);
    lam_vy0 = z0(2);
    lam_y   = z0(3);
    lam_m0  = z0(4);
    t_burn  = z0(5);

    if t_burn <= 0 || t_burn > 2
        res = 1e6 * ones(5,1);
        return;
    end

    pp = p;
    pp.lam_vx0 = lam_vx0;
    pp.lam_vy0 = lam_vy0;
    pp.lam_y   = lam_y;

    ic = [p.x0; p.y0; p.vx0; p.vy0; p.m0; lam_m0];

    try
        [~, Z] = ode45(@(t,z) ode_burn(t, z, pp), [0 t_burn], ic, opts_ode);
        zf = Z(end,:);
    catch
        res = 1e6 * ones(5,1);
        return;
    end

    yc     = zf(2);
    vxc    = zf(3);
    vyc    = zf(4);
    mc     = zf(5);
    lam_mc = zf(6);

    lam_vy_c   = lam_vy0 - lam_y * t_burn;
    lam_v_norm = sqrt(lam_vx0^2 + lam_vy_c^2);

    % Switching function
    S = lam_v_norm / mc - 1 / p.c;

    res = [vxc - 1;
           yc + 0.5 * vyc^2 - p.yf;
           lam_mc - 1;
           lam_vy_c - lam_y * vyc;
           S];
end
