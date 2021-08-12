% look at sensitivity contour from KNM-2 like simulation with flat MTD
% (same qU as KNM-2 MTD)

FakeInitFile = @ref_KNM2_KATRIN_RegMTD;%LinFlatMTD;
range = 40;
freePar = 'E0 Norm Bkg';
nGridSteps = 25;%30;

%% tritium run model
F = RunAnalysis('RunNr',1,...
    'DataType','Fake',...
    'FakeInitFile',FakeInitFile,...
    'fixPar',freePar,...
    'SysBudget',40,...
    'AnaFlag','StackPixel',...
    'chi2','chi2Stat',...
    'FSDFlag','Sibille0p5eV',...
    'ELossFlag','KatrinT2A20',...
    'DopplerEffectFlag','FSD',...
    'RadiativeFlag','ON',...
    'minuitOpt','min ; minos');

F.exclDataStart = F.GetexclDataStart(range);

%% configure Sterile analysis object
SterileArg = {'RunAnaObj',F,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
    'nGridSteps',nGridSteps,...
    'SmartGrid','OFF',...
    'RecomputeFlag','OFF',...
    'SysEffect','KSN2Top4',...
    'RandMC','OFF',...
    'range',range,...
    'LoadGridArg',{'mNu4SqTestGrid','OFF'}};
S = SterileAnalysis(SterileArg{:});

S.RunAnaObj.chi2 = 'chi2Stat';
S.RunAnaObj.NonPoissonScaleFactor = 1;
S.LoadGridFile(S.LoadGridArg{:});
S.Interp1Grid('MaxM4Sq',38.^2);
S.ContourPlot; 
sin2T4_contour_stat = S.sin2T4_contour;
mNu4Sq_contour_stat = S.mNu4Sq_contour;

S.RunAnaObj.chi2 = 'chi2CMShape';
S.RunAnaObj.NonPoissonScaleFactor = 1.112;
 S.GridSearch(S.LoadGridArg{:});
 return;

S.LoadGridFile(S.LoadGridArg{:});
S.Interp1Grid('MaxM4Sq',38.^2);
S.ContourPlot('HoldOn','ON','Color',rgb('Orange'),'LineStyle',':');
sin2T4_contour_cm = S.sin2T4_contour;
mNu4Sq_contour_cm = S.mNu4Sq_contour;
%%
mNu4Sq    = linspace(max([min(mNu4Sq_contour_stat),min(mNu4Sq_contour_cm)]),min([max(mNu4Sq_contour_stat),max(mNu4Sq_contour_cm)]),1e3);
sin2T4_cm = interp1(mNu4Sq_contour_cm,sin2T4_contour_cm,mNu4Sq,'spline');
sin2T4_stat = interp1(mNu4Sq_contour_stat,sin2T4_contour_stat,mNu4Sq,'spline');

%% plot
f2 = figure('Units','normalized','Position',[-0.1,0.1,0.6,0.6]);
s1 = subplot(2,3,[1,2,4,5]);
pstat =plot(sin2T4_contour_stat,mNu4Sq_contour_stat,'Color',rgb('Silver'),'LineWidth',3);
hold on;
ptot =plot(sin2T4_contour_cm,mNu4Sq_contour_cm,'-.','Color',rgb('Black'),'LineWidth',3);
leg = legend([pstat,ptot],'Stat. only','Total (stat. and syst.)','Location','northeast'); 
PrettyLegendFormat(leg); legend boxoff
set(gca,'XScale','log');
set(gca,'YScale','log');
PrettyFigureFormat('FontSize',22);
ax1 = gca;
xlabel(sprintf('|{\\itU}_{e4}|^2'));
ylabel(sprintf('{\\itm}_4^2 (eV)'));
xlim([5e-03 0.5]);

sysOnly = sqrt(sin2T4_cm.^2-sin2T4_stat.^2);

s2 = subplot(2,3,[3,6]);
pStat = plot(sysOnly.^2./sin2T4_cm.^2,mNu4Sq,'-k','LineWidth',2);
hold on;
pTot = plot(sin2T4_stat.^2./sin2T4_cm.^2,mNu4Sq,':','Color',rgb('Silver'),'LineWidth',2);
leg = legend([pStat,pTot],'Stat. only','All syst.','Location','north');

PrettyLegendFormat(leg);
PrettyFigureFormat('FontSize',22);
set(gca,'YScale','log');
ax2 = gca;
ax2.YAxisLocation = 'right';
xlabel(sprintf('\\sigma^2 / \\sigma_{total}^2'));
%xlim([0 1]);
linkaxes([s1,s2],'y');
ax1.Position(3) = 0.55;
ax1.Position(2) = 0.15;
ax2.Position(2) = ax1.Position(2);
ax1.Position(1) = 0.1;
ax2.Position(1) = 0.66;