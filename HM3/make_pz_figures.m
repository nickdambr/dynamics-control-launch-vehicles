%% HM3 - Pole-zero / root-locus companion figures
%  Complementary s-plane views for the classical (Nichols) design:
%    (1) Intro - root locus of the loop gain: conditional stability, the two
%        gain-margin boundaries of the aero-unstable airframe driven through
%        the TVC actuator.
%    (2) Task 1 - closed-loop pole placement: the +sqrt(A6) open-loop pole
%        pulled into the LHP by the rigid PD.
%    (3) Task 2 - bending poles with/without the notch: the lightly damped
%        bending pair (zeta=0.005) driven unstable by the loop and recovered
%        by the minimum-phase notch.
%    (4) Task 3 - closed-loop pole migration over the +/-30% uncertainty box
%        (Nominal + vertices V1..V4): all poles stay in the LHP -> robust.
%
%  Reuses the same plant/controller/loop builders as main_task1/2/3.m, so the
%  s-plane picture is consistent with the Nichols analysis. Writes PNGs to
%  figures/ (owned by this script; shared by both report masters).

clear; close all; clc;
warning('off','Control:analysis:MarginUnstable');

here    = fileparts(mfilename('fullpath'));
fig_dir = fullfile(here,'figures');
if ~exist(fig_dir,'dir'); mkdir(fig_dir); end

% --- house style ---------------------------------------------------------
cUnst = [0.80 0.16 0.16];      % unstable / open-loop / no-notch reference
cStab = [0.00 0.40 0.70];      % stable / closed-loop / with-notch
cLoc  = [0.36 0.36 0.55];      % root-locus branches (single colour)
cNom  = [0.93 0.60 0.13];      % nominal operating point (k=1)
cRHP  = [0.80 0.16 0.16];      % right-half-plane shading (unstable region)
cCorn = [0.00 0.00 0.00;       % Nominal
         0.00 0.45 0.74;       % V1
         0.85 0.33 0.10;       % V2
         0.47 0.67 0.19;       % V3
         0.49 0.18 0.56];      % V4

%% Shared design point (nominal, max-qbar): the report's controller
p0     = load_hw3_params();
Grigid = build_plant_rigid(p0);
K      = design_controller(Grigid, [], 'verbose', false);   % Task 1/2 PD
zN = 0.002; zD = 0.7; sgn = +1;                              % deep min-phase notch

%% (1) Root loci of the loops actually used -- built-in Control System Toolbox
%  rlocus (standard branch rendering, direction and pole/zero conventions),
%  then lightly styled (light theme, RHP shading, nominal k=1 poles). Two loops:
%    - rigid Task-1 loop (ideal actuator): the +sqrt(A6) aero pole is pulled
%      into the LHP -> conditional stability (a MINIMUM loop gain is required);
%    - full Task-2 loop (flex + TVC + delay + notch): zoomed near the origin to
%      show the bending/notch branches (its far branches race toward the RHP
%      zeros of the Pade delay approximation, left off-view).
Wtvc   = build_tvc(p0,3);
Gfull  = build_plant_full(p0,'ins');
notch  = build_notch_filter(p0.wBM, zN, zD, sgn);
[Lrig,~] = assemble_loop(Grigid, K, []);               % Task 1 loop (ideal actuator)
[Lful,~] = assemble_loop(Gfull,  K, Wtvc*notch);       % Task 2/3 loop (TVC+delay+notch)

rl_fig(Lrig, 'Root locus of the rigid Task-1 loop', ...
       'intro_rootlocus', fig_dir, [-10 2.5], [-4 4]);   % wide/tall enough to show the arc pair (|Im|~3.6) break in at s~-7.6 and one branch return to the -4.03 zero
rl_fig(Lful, 'Root locus of the full Task-2 loop (zoom near the origin)', ...
       'task2_rootlocus', fig_dir, [-30 15], [-40 40]);

%% (2) Task 1 - closed-loop pole placement (rigid PD, ideal actuator)
[~,T1] = assemble_loop(Grigid, K, []);
pOL = pole(Grigid);  pCL = pole(T1);

f2 = figure('Name','task1_polemap','Color','w','Position',[100 100 720 560]);
ax = axes(f2); hold(ax,'on');
xl = [-6 2];  yl = [-0.28 0.28];
shade_rhp(ax, xl, yl, cRHP);
hOL = plot(ax, real(pOL),imag(pOL),'x','Color',cUnst,'MarkerSize',13,'LineWidth',2.2);
hCL = plot(ax, real(pCL),imag(pCL),'o','MarkerFaceColor',cStab,'MarkerEdgeColor',cStab*0.55,'MarkerSize',9,'LineWidth',1.0);
finish_ax(ax, xl, yl, 'Task 1: closed-loop pole placement (rigid PD)');
legend([hOL hCL], {'open-loop plant poles','closed-loop poles'}, ...
        'Location','northwest','FontSize',10);

%% (3) Task 2 - bending poles, with vs without the notch
[~,Tno] = assemble_loop(Gfull, K, Wtvc);          % Gfull, notch defined in (1)
[~,Tnf] = assemble_loop(Gfull, K, Wtvc*notch);
pno = pole(Tno);  pnf = pole(Tnf);

f3 = figure('Name','task2_notch_poles','Color','w','Position',[100 100 720 600]);
ax = axes(f3); hold(ax,'on');
xl = [-30 25];  yl = [-45 45];
shade_rhp(ax, xl, yl, cRHP);
yline(ax, p0.wBM,'k:','\omega_{BM}','LabelHorizontalAlignment','left', ...
      'FontSize',11,'HandleVisibility','off');
hNo = plot(ax, real(pno),imag(pno),'x','Color',cUnst,'MarkerSize',12,'LineWidth',2.0);
hWi = plot(ax, real(pnf),imag(pnf),'o','MarkerFaceColor',cStab,'MarkerEdgeColor',cStab*0.55,'MarkerSize',9,'LineWidth',1.0);
finish_ax(ax, xl, yl, 'Task 2: bending poles without / with the notch');
legend([hNo hWi], {'no notch (unstable)','with notch (stable)'}, ...
        'Location','northwest','FontSize',10);

%% (4) Task 3 - closed-loop pole migration over the +/-30% box
cases = {'Nominal',1.00,1.00; 'V1',0.70,0.70; 'V2',0.70,1.30; ...
         'V3',1.30,0.70; 'V4',1.30,1.30};
nC = size(cases,1);

f4 = figure('Name','task3_pole_migration','Color','w','Position',[100 100 720 620]);
ax = axes(f4); hold(ax,'on');
xl = [-30 8];  yl = [-35 35];
shade_rhp(ax, xl, yl, cRHP);
hLeg = gobjects(nC,1);
for i = 1:nC
    p  = load_hw3_params('mu_alpha_scale',cases{i,2},'mu_c_scale',cases{i,3});
    Gf = build_plant_full(p,'ins');
    Wf = build_tvc(p,3) * build_notch_filter(p0.wBM, zN, zD, sgn);   % notch frozen @ nominal
    [~,T] = assemble_loop(Gf, K, Wf);
    pp = pole(T);
    hLeg(i) = plot(ax, real(pp),imag(pp),'x','Color',cCorn(i,:),'MarkerSize',10,'LineWidth',1.7);
end
finish_ax(ax, xl, yl, 'Task 3: closed-loop pole migration over the \pm30% box');
legend(hLeg, cases(:,1), 'Location','northwest','FontSize',10);

%% Export
for f = [f2 f3 f4]
    try
        theme(f,'light');
    catch
        f.Color = 'w';
    end
    exportgraphics(f, fullfile(fig_dir,[get(f,'Name') '.png']),'Resolution',220);
end
fprintf('Pole-zero / root-locus figures written to %s\n', fig_dir);

%% ------------------------------------------------------------------ helpers
function shade_rhp(ax, xl, yl, c)
% Light shading of the right-half-plane (unstable region), behind the data.
if xl(2) > 0
    patch(ax, [0 xl(2) xl(2) 0], [yl(1) yl(1) yl(2) yl(2)], c, ...
          'EdgeColor','none','FaceAlpha',0.08,'HandleVisibility','off');
end
end

function finish_ax(ax, xl, yl, ttl)
% Consistent publication styling: fonts, grid, axes, imaginary axis, labels.
set(ax,'FontName','Helvetica','FontSize',12,'LineWidth',0.9,'Layer','top','Box','on');
xlim(ax,xl); ylim(ax,yl);
grid(ax,'on'); set(ax,'GridColor',[0.80 0.80 0.80],'GridAlpha',1.0);
xline(ax,0,'Color',[0.30 0.30 0.30],'LineWidth',1.2,'HandleVisibility','off');
yline(ax,0,'Color',[0.86 0.86 0.86],'LineWidth',0.8,'HandleVisibility','off');
xlabel(ax,'Real part  (s^{-1})');
ylabel(ax,'Imaginary part  (s^{-1})');
title(ax, ttl, 'FontWeight','bold','FontSize',13);
end

function rl_fig(L, ttl, fname, fig_dir, xl, yl)
% Native Control System Toolbox root locus: rlocusplot(L), with only its own
% settings changed via setoptions (axis limits, title, fonts). No shading, no
% overlays -- exactly the rlocus rendering, just zoomed and titled.
f = figure('Color','w','Position',[100 100 740 560]);
h = rlocusplot(L);
setoptions(h, 'XLim',{xl}, 'YLim',{yl}, 'XLimMode','manual','YLimMode','manual');
opts = getoptions(h);
opts.Title.String  = ttl;   opts.Title.FontSize = 13;
opts.XLabel.FontSize = 12;  opts.YLabel.FontSize = 12;  opts.TickLabel.FontSize = 11;
setoptions(h, opts);
try
    theme(f,'light');
catch
    f.Color = 'w';
end
exportgraphics(f, fullfile(fig_dir,[fname '.png']),'Resolution',220);
end