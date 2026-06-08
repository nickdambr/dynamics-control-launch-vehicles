%% TASK 4 - Optimal staging (two-stage vehicle, no vertical ascent)
%  Two burn arcs with staging (jettison stage 1 structure at ts)
%  Search for optimal staging time to maximize payload

clear; close all; clc;

%% Parameters
c   = 0.6;
eta = 0.1;
yf  = 0.04;
Q   = 2;         % same Q for both stages
T   = c * Q;

opts_ode = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
opts_fs  = optimoptions('fsolve', 'Display', 'off', ...
    'MaxIterations', 500, 'MaxFunctionEvaluations', 5000, ...
    'FunctionTolerance', 1e-10, 'StepTolerance', 1e-10);

%% First solve single-stage (Task 1) as reference
fprintf('=== Single-stage reference solution ===\n');
p1.c = c; p1.Q = Q; p1.T = T; p1.yf = yf;

z_guess_1 = [0.22; 1.3; 4.5; -0.5; 0.30];
[z_ref, ~, ef] = fsolve(@(z) shooting_single(z, p1, opts_ode), z_guess_1, opts_fs);
if ef > 0
    pp = p1;
    pp.lam_vx0 = z_ref(1); pp.lam_vy0 = z_ref(2); pp.lam_y = z_ref(3);
    ic = [0;0;0;0;1;z_ref(4)];
    [~,Z] = ode45(@(t,z) ode_burn(t,z,pp), [0 z_ref(5)], ic, opts_ode);
    mf_single = Z(end,5);
    tf_single = z_ref(5);
    payload_single = mf_single*(1+eta) - eta;
    fprintf('  tf = %.6f,  mf = %.6f,  payload = %.6f\n', tf_single, mf_single, payload_single);
else
    error('Single-stage reference solution did not converge.');
end

%% Sweep staging time to find optimal two-stage solution
fprintf('\n=== Two-stage optimization: sweep staging time ===\n');

% Staging time must be: 0 < ts < tf
% ms1 = eta * Q * ts (structural mass of stage 1)
% After staging: m_plus = m_minus - ms1 = (1 - Q*ts) - eta*Q*ts = 1 - Q*ts*(1+eta)
% This must be > 0: ts < 1/(Q*(1+eta))
ts_max = min(tf_single * 0.95, 1/(Q*(1+eta)) - 0.01);
ts_vec = linspace(0.01, ts_max, 50);

mf_two   = nan(size(ts_vec));
tf_two   = nan(size(ts_vec));
pay_two  = nan(size(ts_vec));
sol_two  = cell(size(ts_vec));

% Use single-stage solution as starting guess
z_prev = [z_ref(1); z_ref(2); z_ref(3); z_ref(4); z_ref(5); z_ref(5)*0.5];
% z0 = [lam_vx0, lam_vy0, lam_y, lam_m0, tf, ts]

% Start from middle of range and sweep outward
[~, idx_mid] = min(abs(ts_vec - tf_single*0.4));

for pass = 1:2
    if pass == 1
        range = idx_mid:length(ts_vec);
    else
        range = idx_mid-1:-1:1;
    end

    z_prev_loc = [];

    for ii = range
        ts = ts_vec(ii);
        p4.c = c; p4.Q = Q; p4.T = T; p4.yf = yf; p4.eta = eta; p4.ts = ts;

        % For fixed ts, solve BVP with unknowns [lam_vx0, lam_vy0, lam_y, lam_m0, tf]
        if ~isempty(z_prev_loc)
            z_guess_4 = z_prev_loc;
        else
            z_guess_4 = z_ref;  % use single-stage solution
        end

        [z_sol, ~, ef] = fsolve(@(z) shooting_twostage(z, p4, opts_ode), z_guess_4, opts_fs);

        if ef > 0
            % Verify and extract solution
            pp = p4;
            pp.lam_vx0 = z_sol(1); pp.lam_vy0 = z_sol(2); pp.lam_y = z_sol(3);
            tf_val = z_sol(5);

            % Integrate stage 1
            ic = [0;0;0;0;1;z_sol(4)];
            [~,Z1] = ode45(@(t,z) ode_burn(t,z,pp), [0 ts], ic, opts_ode);
            z_s = Z1(end,:);

            % Staging
            ms1 = eta * Q * ts;
            z_s(5) = z_s(5) - ms1;  % mass drop

            % Integrate stage 2
            [~,Z2] = ode45(@(t,z) ode_burn(t,z,pp), [ts tf_val], z_s', opts_ode);

            mf_two(ii) = Z2(end,5);
            tf_two(ii) = tf_val;

            % Payload = mf - ms2,  where ms2 = eta * Q * (tf - ts)
            mp2 = Q * (tf_val - ts);
            ms2 = eta * mp2;
            pay_two(ii) = mf_two(ii) - ms2;

            sol_two{ii} = z_sol;
            z_prev_loc = z_sol;
        end
    end
end

%% Find optimal staging
valid = ~isnan(pay_two);
if any(valid)
    [pay_opt, idx_opt] = max(pay_two(valid));
    idx_valid = find(valid);
    idx_opt = idx_valid(idx_opt);
    ts_opt = ts_vec(idx_opt);
    mf_opt = mf_two(idx_opt);
    tf_opt = tf_two(idx_opt);

    fprintf('\n=== Optimal staging results ===\n');
    fprintf('  Staging time:    ts = %.6f\n', ts_opt);
    fprintf('  Total burn time: tf = %.6f\n', tf_opt);
    fprintf('  Final mass:      mf = %.6f\n', mf_opt);
    fprintf('  Payload (2-st):  %.6f\n', pay_opt);
    fprintf('  Payload (1-st):  %.6f\n', payload_single);
    fprintf('  Payload gain:    %.6f (%.2f%%)\n', ...
        pay_opt - payload_single, (pay_opt/payload_single - 1)*100);

    %% Plot: Payload vs staging time
    figure('Name','Task 4 - Payload vs staging time');
    plot(ts_vec(valid), pay_two(valid), 'b-o', 'MarkerSize', 3);
    hold on;
    yline(payload_single, 'r--', 'Single stage', 'LineWidth', 1.5);
    plot(ts_opt, pay_opt, 'rp', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
    xlabel('Staging time t_s (nondim)');
    ylabel('Payload mass (nondim)');
    title('Task 4: Payload vs staging time');
    grid on; legend('Two-stage', 'Single-stage', 'Optimal', 'Location','best');

    %% Plot: Final mass vs staging time
    figure('Name','Task 4 - Final mass vs staging time');
    plot(ts_vec(valid), mf_two(valid), 'b-o', 'MarkerSize', 3);
    hold on;
    yline(mf_single, 'r--', 'Single stage m_f', 'LineWidth', 1.5);
    xlabel('Staging time t_s (nondim)');
    ylabel('Final mass m_f (nondim)');
    title('Task 4: Final mass vs staging time');
    grid on;

    %% Plot optimal trajectory
    z_opt = sol_two{idx_opt};
    pp = struct('c', c, 'Q', Q, 'T', T, 'yf', yf, 'eta', eta, 'ts', ts_opt);
    pp.lam_vx0 = z_opt(1); pp.lam_vy0 = z_opt(2); pp.lam_y = z_opt(3);

    ic = [0;0;0;0;1;z_opt(4)];
    [T1, Z1] = ode45(@(t,z) ode_burn(t,z,pp), linspace(0, ts_opt, 300), ic, opts_ode);

    z_s = Z1(end,:);
    ms1 = eta * Q * ts_opt;
    z_s(5) = z_s(5) - ms1;
    [T2, Z2] = ode45(@(t,z) ode_burn(t,z,pp), linspace(ts_opt, tf_opt, 300), z_s', opts_ode);

    figure('Name','Task 4 - Optimal trajectory');
    plot(Z1(:,1), Z1(:,2), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Stage 1');
    hold on;
    plot(Z2(:,1), Z2(:,2), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Stage 2');
    plot(Z2(1,1), Z2(1,2), 'ko', 'MarkerSize', 8, 'MarkerFaceColor','k', ...
        'DisplayName', 'Staging point');
    xlabel('x (nondim)'); ylabel('y (nondim)');
    title('Task 4: Optimal two-stage trajectory');
    legend('Location','best'); grid on;

    % Mass profile
    figure('Name','Task 4 - Mass profile');
    plot(T1, Z1(:,5), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Stage 1');
    hold on;
    plot([ts_opt ts_opt], [Z1(end,5) z_s(5)], 'k--', 'LineWidth', 1.5, ...
        'DisplayName', 'Staging (jettison)');
    plot(T2, Z2(:,5), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Stage 2');
    xlabel('Time (nondim)'); ylabel('Mass (nondim)');
    title('Task 4: Mass profile');
    legend('Location','best'); grid on;
else
    fprintf('ERROR: No converged two-stage solutions found.\n');
end

fprintf('\nTask 4 complete.\n');

%% ===================== EXPORT FIGURES =====================
fig_dir = fullfile(fileparts(mfilename('fullpath')), 'figures');
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end
slugify = @(s) lower(regexprep(s, '[^a-zA-Z0-9]+', '_'));
fig_handles = findobj(groot, 'Type', 'figure');
for kk = 1:numel(fig_handles)
    nm = get(fig_handles(kk), 'Name');
    if isempty(nm); nm = sprintf('fig%d', kk); end
    try
        theme(fig_handles(kk), 'light');    % force light theme (ignore desktop dark mode)
        drawnow;
    catch
        set(fig_handles(kk), 'Color', 'w'); % fallback for pre-R2025a MATLAB
    end
    exportgraphics(fig_handles(kk), ...
        fullfile(fig_dir, ['task4_' slugify(nm) '.png']), 'Resolution', 200);
end

%% ===================== LOCAL FUNCTIONS =====================

function res = shooting_single(z0, p, opts_ode)
% Single-stage shooting (same as Task 1)
    lam_vx0 = z0(1); lam_vy0 = z0(2); lam_y = z0(3);
    lam_m0 = z0(4);  tf = z0(5);

    if tf <= 0 || tf > 2
        res = 1e6*ones(5,1); return;
    end

    pp = p;
    pp.lam_vx0 = lam_vx0; pp.lam_vy0 = lam_vy0; pp.lam_y = lam_y;

    ic = [0;0;0;0;1;lam_m0];
    try
        [~,Z] = ode45(@(t,z) ode_burn(t,z,pp), [0 tf], ic, opts_ode);
        zf = Z(end,:);
    catch
        res = 1e6*ones(5,1); return;
    end

    lam_vy_f = lam_vy0 - lam_y * tf;
    lam_v_norm = sqrt(lam_vx0^2 + lam_vy_f^2);
    H_f = lam_y*zf(4) + (p.T/zf(5))*lam_v_norm - lam_vy_f - p.Q*zf(6);

    res = [zf(2)-p.yf; zf(3)-1; zf(4); zf(6)-1; H_f];
end

function res = shooting_twostage(z0, p, opts_ode)
% Two-stage shooting with fixed staging time p.ts
%   z0 = [lam_vx0; lam_vy0; lam_y; lam_m0; tf]

    lam_vx0 = z0(1); lam_vy0 = z0(2); lam_y = z0(3);
    lam_m0 = z0(4);  tf = z0(5);
    ts = p.ts;

    if tf <= ts || tf > 2 || ts <= 0
        res = 1e6*ones(5,1); return;
    end

    pp = p;
    pp.lam_vx0 = lam_vx0; pp.lam_vy0 = lam_vy0; pp.lam_y = lam_y;

    % Stage 1: integrate from 0 to ts
    ic = [0;0;0;0;1;lam_m0];
    try
        [~,Z1] = ode45(@(t,z) ode_burn(t,z,pp), [0 ts], ic, opts_ode);
        z_s = Z1(end,:)';
    catch
        res = 1e6*ones(5,1); return;
    end

    % Staging: drop structure of stage 1
    ms1 = p.eta * p.Q * ts;
    z_s(5) = z_s(5) - ms1;
    % lam_m is continuous across staging
    if z_s(5) <= 0
        res = 1e6*ones(5,1); return;
    end

    % Stage 2: integrate from ts to tf
    try
        [~,Z2] = ode45(@(t,z) ode_burn(t,z,pp), [ts tf], z_s, opts_ode);
        zf = Z2(end,:);
    catch
        res = 1e6*ones(5,1); return;
    end

    lam_vy_f = lam_vy0 - lam_y * tf;
    lam_v_norm = sqrt(lam_vx0^2 + lam_vy_f^2);
    H_f = lam_y*zf(4) + (p.T/zf(5))*lam_v_norm - lam_vy_f - p.Q*zf(6);

    res = [zf(2)-p.yf; zf(3)-1; zf(4); zf(6)-1; H_f];
end
