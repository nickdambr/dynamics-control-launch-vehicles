%% HM3 - Monte-Carlo stability-margin robustness (probabilistic Task 3)
%  Probabilistic counterpart of main_task3.m. The Task-2 controller (rigid PD
%  + bending notch fixed at the NOMINAL bending frequency) is held fixed, and
%  the physical uncertainties are now treated as RANDOM VARIABLES rather than
%  the four deterministic +/-30 % box vertices. For each Monte-Carlo draw the
%  open-loop transfer L is reassembled and its stability margins are recomputed,
%  yielding distributions of gain/phase/delay margin, a probability of stability
%  P(stable), a Nichols "cloud", and the parameter sensitivities that drive the
%  worst cases.
%
%  Uncertain parameters (multiplicative factors on the nominal Table-1 values):
%     mu_alpha = A6   aerodynamic moment        Gaussian, 3-sigma = +/-30 % (Task-3 box)
%     mu_c     = K1   control effectiveness     Gaussian, 3-sigma = +/-30 % (Task-3 box)
%     wBM             first bending frequency   Gaussian, 3-sigma = +/-6 %  (notch detuning)
%     zBM             bending damping           lognormal, factor ~ [0.45, 2.2] at 2-sigma
%     tau             TVC transport delay       uniform,  +/-25 %
%
%  The notch is NOT retuned when wBM disperses: this is exactly the realistic
%  robustness question (the deep notch needs near-exact wBM knowledge, see the
%  Task-2 detuning note in the README), so wBM is the dominant driver.
%
%  Dependencies: Control System Toolbox only. parfor degrades to a serial loop
%  without the Parallel Computing Toolbox; percentiles use a base-MATLAB helper
%  (no Statistics Toolbox needed).
%
%  Reference: Homework 3 (Zavoli, v1.2, May 2026), Task 3 (extension).

clear; close all; clc;
warning('off','Control:analysis:MarginUnstable');

%% ----------------------------------------------------------- Configuration
N      = 1500;                 % number of Monte-Carlo samples (bump to 5e3+ if PCT available)
Nsub   = 150;                  % loop transfers overlaid in the Nichols cloud
wgrid  = logspace(-2, 2, 600); % rad/s, frequency grid for the cloud
seed   = 2026;                 % RNG seed (reproducible)

% Uncertainty specification (multiplicative factor on each nominal coefficient)
unc.mu_alpha = struct('dist','gauss',  'sigma',0.10, 'trunc',0.30);  % A6
unc.mu_c     = struct('dist','gauss',  'sigma',0.10, 'trunc',0.30);  % K1
unc.wBM      = struct('dist','gauss',  'sigma',0.02, 'trunc',0.06);  % bending freq
unc.zBM      = struct('dist','lognorm','sigma',0.40);                % bending damping
unc.tau      = struct('dist','uniform','half', 0.25);                % TVC delay

%% --------------------------- Fixed controller (Task-2 design) + nominal model
p0     = load_hw3_params();
Grigid = build_plant_rigid(p0);
K      = design_controller(Grigid, [], 'verbose', false);
notch  = struct('wx',p0.wBM,'zN',0.002,'zD',0.7,'sgn',+1);   % fixed @ nominal wBM
Wnotch = build_notch_filter(notch.wx,notch.zN,notch.zD,notch.sgn);
w      = load_wind_profile(p0);                              % same gust for all draws

% Nominal loop (sanity reference + cloud baseline)
Gn      = build_plant_full(p0,'ins');
Ln      = assemble_loop(Gn, K, build_tvc(p0,3)*Wnotch);
[~, mn] = nichols_branch(Ln, wgrid, []);   % nominal phase = cloud alignment reference

fprintf('Fixed controller: Kp_th=%.3f Kd_th=%.3f | notch wx=%.1f zN=%.3f zD=%.2f\n', ...
        K.Kp_th, K.Kd_th, notch.wx, notch.zN, notch.zD);
fprintf('Monte-Carlo: N=%d samples, seed=%d\n', N, seed);

%% --------------------------------------------------- Pre-sample uncertainties
%  Sampling is done OUTSIDE the parfor (deterministic given the seed, and keeps
%  the parallel loop body RNG-free).
rng(seed);
fa = sample_factor(unc.mu_alpha, N);   % mu_alpha = A6 factor
fc = sample_factor(unc.mu_c,     N);   % mu_c     = K1 factor
fw = sample_factor(unc.wBM,      N);   % bending-frequency factor
fz = sample_factor(unc.zBM,      N);   % bending-damping factor
ft = sample_factor(unc.tau,      N);   % TVC-delay factor

%% ----------------------------------------------------- Monte-Carlo propagation
rigidGM = nan(N,1);   % low-frequency aerodynamic gain margin [dB]
minGM   = nan(N,1);   % worst gain margin over all crossovers [dB]
PM      = nan(N,1);   % phase margin [deg]
DM      = nan(N,1);   % delay margin [ms]
peakTh  = nan(N,1);   % peak pitch attitude under the gust [deg]
peakZ   = nan(N,1);   % peak lateral drift under the gust [m]
stab    = false(N,1); % closed-loop stability flag

t_mc = tic;
parfor i = 1:N
    % perturbed parameters (start from nominal, scale the uncertain fields)
    p = p0;
    p.A6  = p0.A6  * fa(i);
    p.K1  = p0.K1  * fc(i);
    p.wBM = p0.wBM * fw(i);
    p.zBM = p0.zBM * fz(i);
    p.tau = p0.tau * ft(i);

    % reassemble the loop with the FIXED controller + FIXED (nominal) notch
    Gf      = build_plant_full(p,'ins');
    Wf      = build_tvc(p,3) * Wnotch;
    [L, T]  = assemble_loop(Gf, K, Wf);

    % --- frequency-domain margins ---
    am  = allmargin(L);
    gmv = 20*log10(am.GainMargin);          % all gain margins [dB]
    gf  = am.GMFrequency;

    if isempty(gmv), minGM(i) = Inf;        % no gain crossover => infinite GM
    else,            minGM(i) = min(abs(gmv)); end

    idx = find(gf>0.2 & gf<1, 1);           % low-freq aerodynamic crossover
    if ~isempty(idx), rigidGM(i) = abs(gmv(idx)); end

    if isempty(am.DelayMargin), DM(i) = Inf;
    else,                       DM(i) = min(am.DelayMargin)*1000; end   % ms

    isStab  = isstable(T);
    stab(i) = isStab;

    [~,Pm] = margin(L);
    if isnan(Pm)                            % no phase crossover
        PM(i) = 180*double(isStab);         % stable -> effectively unbounded; cap at 180
    else
        PM(i) = abs(Pm);
    end

    % --- time-domain dispersion (gust response) ---
    r         = simulate_gust_response(T, w);
    peakTh(i) = r.peak_theta*180/pi;
    peakZ(i)  = r.peak_z;
end
fprintf('Monte-Carlo done in %.1f s.\n', toc(t_mc));

%% ----------------------------------------------------------------- Statistics
Pstab  = mean(stab);
Pgm3   = mean(minGM >= 3);         % >= 3 dB design-margin guideline
Ppm30  = mean(PM   >= 30);         % assignment phase-margin magnitude target
qs     = [5 50 95];                % reported percentiles

[~,iWorst] = min(minGM);           % tightest-margin draw

fprintf('\n================  Monte-Carlo robustness summary  ================\n');
fprintf('P(closed-loop stable)        = %.1f %%  (%d/%d)\n', 100*Pstab, sum(stab), N);
fprintf('P(min|GM| >= 3 dB)           = %.1f %%\n', 100*Pgm3);
fprintf('P(|PM| >= 30 deg)            = %.1f %%\n', 100*Ppm30);
fprintf('(|GM|,|PM| are magnitudes: the airframe is open-loop unstable so the\n');
fprintf(' loop is conditionally stable -- the binding stability indicator is the\n');
fprintf(' isstable() flag, not the sign of the gain margin.)\n');
fprintf('\n%-12s %8s %8s %8s\n','metric','p05','p50','p95');
print_pct('min|GM| [dB]', minGM,  qs);
print_pct('rigid|GM|[dB]',rigidGM,qs);
print_pct('|PM| [deg]',   PM,     qs);
print_pct('DM [ms]',      DM,     qs);
print_pct('peak th [deg]',peakTh, qs);
print_pct('peak z [m]',   peakZ,  qs);
fprintf('\nWorst draw (min|GM|=%.2f dB, stable=%d): mu_a=%.2f mu_c=%.2f ', ...
        minGM(iWorst), stab(iWorst), fa(iWorst), fc(iWorst));
fprintf('dwBM=%+.1f%% zBM x%.2f tau x%.2f\n', 100*(fw(iWorst)-1), fz(iWorst), ft(iWorst));
fprintf('===================================================================\n');

%% ----------------------------------------------------------------- Figure 1: histograms
f1 = figure('Name','margins_hist','Color','w','Position',[80 80 1000 560]);
tl = tiledlayout(f1,2,3,'TileSpacing','compact','Padding','compact');
title(tl, sprintf('HM3 Monte-Carlo (N=%d): margin distributions  —  P(stable)=%.1f%%', N, 100*Pstab));

hist_metric(minGM,   'min |GM| [dB]',   [], 'b');          % magnitude (conditional stability)
hist_metric(rigidGM, 'rigid |GM| [dB]', [6 NaN], 'r');     % assignment 6 dB target
hist_metric(PM,      '|PM| [deg]',      [30 NaN],'r');     % assignment 30 deg target
hist_metric(DM,      'delay margin [ms]', [], 'b');
hist_metric(peakTh,  'peak \theta [deg]', [], 'b');
hist_metric(peakZ,   'peak z [m]',        [], 'b');

%% ----------------------------------------------------------------- Figure 2: Nichols cloud
f2 = figure('Name','nichols_cloud','Color','w','Position',[120 120 680 580]);
ax = axes(f2); hold(ax,'on'); grid(ax,'on'); box(ax,'on');

% nominal branch + a single global shift so its 0-dB crossover sits near -180 deg
[gn, phn] = nichols_branch(Ln, wgrid, mn);
sc = find(gn(1:end-1).*gn(2:end) <= 0, 1, 'last');         % last 0-dB gain crossover
if isempty(sc), shift = 360*round(median(phn)/360);
else,           shift = 360*round((phn(sc)+180)/360); end

sub = round(linspace(1, N, min(Nsub,N)));     % representative subset
for k = sub
    p = p0;
    p.A6=p0.A6*fa(k); p.K1=p0.K1*fc(k); p.wBM=p0.wBM*fw(k);
    p.zBM=p0.zBM*fz(k); p.tau=p0.tau*ft(k);
    Lk = assemble_loop(build_plant_full(p,'ins'), K, build_tvc(p,3)*Wnotch);
    [g, ph] = nichols_branch(Lk, wgrid, mn);   % aligned to nominal branch
    if stab(k), col = [0.55 0.55 0.55 0.18]; else, col = [0.85 0.20 0.20 0.35]; end
    plot(ax, ph - shift, g, 'Color', col, 'LineWidth', 0.5);
end
plot(ax, phn - shift, gn, 'b', 'LineWidth', 2);            % nominal
yline(ax, 0, 'k:');                                        % 0 dB
xl = xlim(ax);
for cp = -900:360:540                                       % critical points (-180 + 360k)
    if cp>=xl(1) && cp<=xl(2), plot(ax, cp, 0, 'r+', 'MarkerSize',10,'LineWidth',1.4); end
end
xlabel(ax,'open-loop phase [deg]'); ylabel(ax,'open-loop gain [dB]');
title(ax, sprintf('Nichols cloud (%d draws) — gray: stable, red: unstable', numel(sub)));

%% ----------------------------------------------------------------- Figure 3: sensitivity scatter
f3 = figure('Name','sensitivity','Color','w','Position',[160 160 980 400]);
tl3 = tiledlayout(f3,1,2,'TileSpacing','compact','Padding','compact');
title(tl3,'Monte-Carlo sensitivities of the worst gain margin');

% (a) (mu_alpha, mu_c) plane: stable vs unstable draws, with the Task-3 box + vertices.
%     Task 3 found all four box vertices stable in (mu_alpha,mu_c) alone; once the
%     bending/delay uncertainties are dispersed too, instability appears inside the box.
ax1 = nexttile; hold(ax1,'on'); box(ax1,'on');
plot(ax1, fa(stab),  fc(stab),  '.', 'Color',[0.20 0.55 0.45],'MarkerSize',7);  % stable
plot(ax1, fa(~stab), fc(~stab), 'rx','MarkerSize',6,'LineWidth',1.0);            % unstable
rectangle(ax1,'Position',[0.7 0.7 0.6 0.6],'EdgeColor','k','LineStyle','--','LineWidth',1.1);
Vx=[0.7 0.7 1.3 1.3]; Vy=[0.7 1.3 0.7 1.3]; lab={'V1','V2','V3','V4'};
plot(ax1, Vx, Vy, 'ks','MarkerFaceColor','k','MarkerSize',6);
text(ax1, Vx+0.012, Vy, lab, 'FontSize',8);
xlabel(ax1,'\mu_\alpha factor'); ylabel(ax1,'\mu_c factor');
title(ax1,'(a) stability over the Task-3 \pm30% box');
legend(ax1,{'stable','unstable'},'Location','best');

% (b) bending-frequency detuning is the dominant driver
ax2 = nexttile; hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
fin = isfinite(minGM);
gs = fin & stab; gu = fin & ~stab;
plot(ax2, 100*(fw(gs)-1), minGM(gs), '.', 'Color',[0.4 0.4 0.4],'MarkerSize',7);
plot(ax2, 100*(fw(gu)-1), minGM(gu), 'rx','MarkerSize',7,'LineWidth',1.2);
yline(ax2, 0, 'r--','LineWidth',1.1);
xlabel(ax2,'bending-freq dispersion \Delta\omega_{BM} [%]'); ylabel(ax2,'min |GM| [dB]');
title(ax2,'(b) notch detuning vs worst margin');
legend(ax2,{'stable','unstable'},'Location','best');

%% ----------------------------------------------------------------- Export
fig_dir = fullfile(fileparts(mfilename('fullpath')), 'figures');
if ~exist(fig_dir,'dir'); mkdir(fig_dir); end
for f = [f1 f2 f3]
    try
        theme(f,'light');
    catch
        f.Color = 'w';
    end
    exportgraphics(f, fullfile(fig_dir, ['mc_' get(f,'Name') '.png']), 'Resolution', 200);
end

% Save the raw results for the report tables (factors + per-draw metrics)
MC = struct('N',N,'seed',seed,'unc',unc, ...
            'fa',fa,'fc',fc,'fw',fw,'fz',fz,'ft',ft, ...
            'rigidGM',rigidGM,'minGM',minGM,'PM',PM,'DM',DM, ...
            'peakTh',peakTh,'peakZ',peakZ,'stab',stab, ...
            'Pstab',Pstab,'Pgm3',Pgm3,'Ppm30',Ppm30);
save(fullfile(fileparts(mfilename('fullpath')),'mc_results.mat'),'MC');
fprintf('\nFigures written to %s\nResults saved to mc_results.mat\n', fig_dir);

%% ============================================================ local functions
function f = sample_factor(spec, n)
%SAMPLE_FACTOR  Draw n multiplicative uncertainty factors (median/mean ~ 1).
switch lower(spec.dist)
    case 'gauss'
        f = 1 + spec.sigma*randn(n,1);
        if isfield(spec,'trunc') && ~isempty(spec.trunc)   % resample outside the box
            lo = 1-spec.trunc; hi = 1+spec.trunc;
            bad = f<lo | f>hi;
            while any(bad)
                f(bad) = 1 + spec.sigma*randn(nnz(bad),1);
                bad = f<lo | f>hi;
            end
        end
    case 'uniform'
        f = 1 + spec.half*(2*rand(n,1)-1);
    case 'lognorm'
        f = exp(spec.sigma*randn(n,1));
    otherwise
        error('sample_factor:dist','unknown distribution ''%s''.',spec.dist);
end
end

function [g, ph] = nichols_branch(L, w, ref)
%NICHOLS_BRANCH  Aligned Nichols data (gain dB, phase deg) on grid w.
%   When ref (a nominal phase vector) is given, the curve is shifted by a
%   multiple of 360 deg so its low-frequency branch matches the nominal one,
%   which keeps the Monte-Carlo cloud on a single branch.
[mag, phase] = nichols(L, w);
g  = 20*log10(squeeze(mag));
ph = squeeze(phase);
if ~isempty(ref)
    ph = ph - 360*round((ph(1)-ref(1))/360);
end
end

function print_pct(name, x, q)
%PRINT_PCT  Print 5/50/95 percentiles of finite samples (base-MATLAB).
v = local_quantile(x, q);
fprintf('%-12s %8.2f %8.2f %8.2f\n', name, v(1), v(2), v(3));
end

function qv = local_quantile(x, p)
%LOCAL_QUANTILE  Percentiles without the Statistics Toolbox (midpoint rule).
x = sort(x(isfinite(x)));
if isempty(x), qv = nan(size(p)); return; end
n  = numel(x);
pos = min(max(p/100*n + 0.5, 1), n);
qv  = interp1(1:n, x, pos, 'linear');
end

function hist_metric(x, name, refline, refcol)
%HIST_METRIC  One histogram tile with optional reference line + median.
nexttile; hold on; grid on; box on;
v = x(isfinite(x));
histogram(v, 'NumBins', 30, 'FaceColor',[0.30 0.50 0.75], 'EdgeColor','none');
med = median(v,'omitnan');
xline(med, 'k-', sprintf('p50=%.2g',med), 'LabelOrientation','horizontal');
if ~isempty(refline) && ~isnan(refline(1))
    xline(refline(1), [refcol '--'], 'LineWidth',1.2);
end
xlabel(name); ylabel('count');
end
