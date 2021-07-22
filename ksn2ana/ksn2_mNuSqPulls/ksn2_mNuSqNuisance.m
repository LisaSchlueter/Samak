% ksn2 calculate chi2 grid search
% compare m2 free, m2 nuisance parameter
%% settings that might change
chi2 = 'chi2CMShape';
DataType = 'Twin';%'Real';
nGridSteps = 40;
range = 40;
Mode = 'Compute';
%% configure RunAnalysis object
if strcmp(chi2,'chi2Stat')
    NonPoissonScaleFactor = 1;
elseif  strcmp(chi2,'chi2CMShape')
    NonPoissonScaleFactor = 1.112;
end
RunAnaArg = {'RunList','KNM2_Prompt',...
    'DataType',DataType,...
    'fixPar','mNu E0 Norm Bkg',...%free par
    'SysBudget',40,...
    'fitter','minuit',...
    'minuitOpt','min;migrad',...
    'RadiativeFlag','ON',...
    'FSDFlag','KNM2_0p5eV',...
    'ELossFlag','KatrinT2A20',...
    'AnaFlag','StackPixel',...
    'chi2',chi2,...
    'NonPoissonScaleFactor',NonPoissonScaleFactor,...
    'FSD_Sigma',sqrt(0.0124+0.0025),...
    'TwinBias_FSDSigma',sqrt(0.0124+0.0025),...
    'TwinBias_Q',18573.7,...
    'PullFlag',99,...;%99 = no pull
    'BKG_PtSlope',3*1e-06,...
    'TwinBias_BKG_PtSlope',3*1e-06,...
    'DopplerEffectFlag','FSD'};
A = MultiRunAnalysis(RunAnaArg{:});
%% configure Sterile analysis object
SterileArg = {'RunAnaObj',A,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
    'nGridSteps',nGridSteps,...
    'SmartGrid','OFF',...
    'RecomputeFlag','OFF',...
    'SysEffect','all',...
    'RandMC','OFF',...
    'range',range,...
    'LoadGridArg',{'mNu4SqTestGrid',5}};

S = SterileAnalysis(SterileArg{:});
S.InterpMode = 'spline';
%
S.LoadGridFile(S.LoadGridArg{:}); 
S.Interp1Grid('nInter',1e3);
%
PlotPar = S.mNuSq;

if strcmp(DataType,'Real')
    ContourVec = [-5 -1,0,0.3,0.9,1,2,10];
    LabelSpacing = 380;
      BF = 'ON';
else
    ContourVec = [-5 -1,-0.3 0,0.3,1,2,10];
    LabelSpacing = 10000;
    BF = 'OFF';
end
Splines = 'ON';
GetFigure;
%
for i=1:numel(ContourVec)
    [M,p1] = contour3(S.sin2T4,S.mNu4Sq,PlotPar,[ContourVec(i),ContourVec(i)],...
        'Color',rgb('LightGreen'),'ShowText','off','LineWidth',2,'LabelSpacing',LabelSpacing);
    hold on;
   % cl = clabel(M,p1,'FontSize',18,'FontName','Times New Roman');
end
%

[pFree,pFix] = S.PlotmNuSqOverview('PullmNuSq','OFF','SavePlot','OFF','HoldOn','ON','BestFit',BF);
view(2)
grid off
pFix.LineWidth = 2;
pFree.LineWidth = 3;
leg = legend([pFix,pFree,p1],...
    sprintf('i)  Fixed {\\itm}_\\nu^2 = 0 eV^2'),...
    sprintf('ii) Free {\\itm}_\\nu^2 unconstrained'),...
    sprintf('Isoline: {\\itm}_\\nu^2 best fit for ii)'),...
    'EdgeColor',rgb('Silver'),'Location','southwest','box','off');
if strcmp(DataType,'Twin')
    legend boxon
    PrettyLegendFormat(leg,'alpha',0.5);
end
PRLFormat;
title('');
ax = gca;
ax.FontSize =  24;
leg.FontSize = 24;
ax.XLabel.FontSize = 28;
ax.YLabel.FontSize = 28;
ylim([1 40^2]);
xlim([2e-03 0.5]);
set(gca,'YScale','log')
set(gca,'XScale','log')
%% save
name_i = strrep(S.DefPlotName,'_mNuE0BkgNorm','');
plotname = sprintf('%s_mNuSqOverviewmNuSq_%.2gCL.png',name_i,S.ConfLevel);
print(gcf,plotname,'-dpng','-r450');
fprintf('save plot to %s \n',plotname);
ylabel(sprintf('{\\itm}_4^2 (eV^{ 2})'));
export_fig(strrep(plotname,'.png','.pdf'));