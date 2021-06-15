% ksn2 calculate chi2 grid search
% compare m2 free, m2 nuisance parameter
%% settings that might change
chi2 = 'chi2Stat';
DataType = 'Twin';
nGridSteps = 25;
range = 40;
%% configure RunAnalysis object
RunAnaArg = {'RunList','KNM2_Prompt',...
    'DataType',DataType,...
    'fixPar','E0 Norm Bkg',...%free par
    'SysBudget',40,...
    'fitter','minuit',...
    'minuitOpt','min;migrad',...
    'RadiativeFlag','ON',...
    'FSDFlag','KNM2_0p5eV',...
    'ELossFlag','KatrinT2A20',...
    'AnaFlag','StackPixel',...
    'chi2',chi2,...
    'NonPoissonScaleFactor',1,...
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
    'range',range};

%%
S = SterileAnalysis(SterileArg{:});
S.InterpMode = 'spline';
pHandle = cell(5,1);

S.nGridSteps = 30;
S.RunAnaObj.chi2 = 'chi2CMShape';
S.RunAnaObj.NonPoissonScaleFactor = 1.112;
S.LoadGridFile('IgnoreKnm2FSDbinning','ON','mNu4SqTestGrid',5,'ExtmNu4Sq','OFF');
S.Interp1Grid('Maxm4Sq',38.2^2);
pHandle{5} =  S.ContourPlot('HoldOn','OFF','SavePlot','OFF','Color',rgb('SlateGray'),'LineStyle','-.');

S.RunAnaObj.chi2 = 'chi2CMShape';
S.nGridSteps = 30;
S.SysEffect = 'BkgPT';
S.RunAnaObj.NonPoissonScaleFactor = 1;
S.LoadGridFile('mNu4SqTestGrid',5);
S.Interp1Grid('Maxm4Sq',38.2^2);
pHandle{2} =  S.ContourPlot('HoldOn','ON','SavePlot','OFF','Color',rgb('FireBrick'),'LineStyle',':');

 S.RunAnaObj.chi2 = 'chi2Stat';
 S.LoadGridFile('IgnoreKnm2FSDbinning','ON','mNu4SqTestGrid',5,'ExtmNu4Sq','ON');
 S.Interp1Grid('Maxm4Sq',38.2^2);
 pHandle{1} =  S.ContourPlot('HoldOn','ON','SavePlot','OFF','Color',rgb('Black'),'LineStyle','-');
 
 S.nGridSteps = 25;
 S.RunAnaObj.NonPoissonScaleFactor = 1.112;
 S.LoadGridFile('IgnoreKnm2FSDbinning','ON','mNu4SqTestGrid','OFF','ExtmNu4Sq','OFF');
 S.Interp1Grid('Maxm4Sq',38.2^2);
 pHandle{3} =  S.ContourPlot('HoldOn','ON','SavePlot','OFF','Color',rgb('LimeGreen'),'LineStyle','--');
 
 S.nGridSteps = 25;
 S.RunAnaObj.NonPoissonScaleFactor = 1;
 S.RunAnaObj.chi2 = 'chi2CMShape';
 S.SysEffect = 'LongPlasma';
 S.LoadGridFile('IgnoreKnm2FSDbinning','ON','mNu4SqTestGrid','OFF','ExtmNu4Sq','OFF');
 S.Interp1Grid('Maxm4Sq',38.2^2);
 pHandle{4} =  S.ContourPlot('HoldOn','ON','SavePlot','OFF','Color',rgb('DodgerBlue'),'LineStyle','-.');
 
% HoldOn = 'ON';

leg = legend([pHandle{:}],'Stat. only',sprintf('Stat. and {\\itB}_{pt}'),...
    sprintf('Stat. and {\\itB}_{NP}'),sprintf('Stat. and Plasma'),'Total');
PrettyLegendFormat(leg);

xlim([0.009,0.012])
ylim([150,620])

%% save
plotname999 = sprintf('%s_BkgPT.png',extractBefore(S.DefPlotName,'FSD'),A.SysBudget);
print(gcf,plotname999,'-dpng','-r300');
fprintf('save plot to %s \n',plotname999)

%%
S.SysEffect = 'all';
S.nGridSteps = 30;
S.RunAnaObj.chi2 = 'chi2CMShape';
S.RunAnaObj.NonPoissonScaleFactor = 1.112;

S.LoadGridFile('IgnoreKnm2FSDbinning','ON','mNu4SqTestGrid',5,'ExtmNu4Sq','OFF');
S.Interp1Grid('Maxm4Sq',38.2^2);
pNew =  S.ContourPlot('HoldOn','OFF','SavePlot','OFF','Color',rgb('Black'),'LineStyle','-');

S.nGridSteps = 50;
S.LoadGridFile('IgnoreKnm2FSDbinning','ON','ExtmNu4Sq','ON');
S.Interp1Grid('Maxm4Sq',38.2^2);
pOld =  S.ContourPlot('HoldOn','ON','SavePlot','OFF','Color',rgb('SlateGray'),'LineStyle','-.');

leg = legend([pNew,pOld],sprintf('Total with {\\itB}_{pt}'),...
          sprintf('Total without {\\itB}_{pt}'));
PrettyLegendFormat(leg);

xlim([5e-03,0.5])
ylim([1,40^2]);

plotname999 = sprintf('%s_TotalBkgPT.png',extractBefore(S.DefPlotName,'FSD'),A.SysBudget);
print(gcf,plotname999,'-dpng','-r300');
fprintf('save plot to %s \n',plotname999)
%%
xlim([0.0098,0.012])
ylim([180,620])

plotname999 = sprintf('%s_TotalBkgPTZoom.png',extractBefore(S.DefPlotName,'FSD'),A.SysBudget);
print(gcf,plotname999,'-dpng','-r300');
fprintf('save plot to %s \n',plotname999)