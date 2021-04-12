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
pHandle = cell(3,1);

% S.nGridSteps = 50;
% S.RunAnaObj.chi2 = 'chi2CMShape';
% S.RunAnaObj.NonPoissonScaleFactor = 1.112;
% S.LoadGridFile('IgnoreKnm2FSDbinning','ON');
% S.Interp1Grid('Maxm4Sq',38.2^2);
% pHandle{1} =  S.ContourPlot('HoldOn','OFF','SavePlot','OFF','Color',rgb('DodgerBlue'),'LineStyle','-');

S.RunAnaObj.chi2 = 'chi2CMShape';
S.nGridSteps = 30;
S.SysEffect = 'BkgPT';
S.RunAnaObj.NonPoissonScaleFactor = 1;
S.GridSearch('mNu4SqTestGrid',5);

return;
%S.LoadGridFile('CheckSmallerN','ON')
%S.Interp1Grid('Maxm4Sq',36^2);
% pHandle{3} =  S.ContourPlot('HoldOn','ON','SavePlot','OFF','Color',rgb('FireBrick'),'LineStyle',':');
% fname = S.GridFilename;
% 
% S.nGridSteps = 50;
% S.RunAnaObj.chi2 = 'chi2Stat';
% S.LoadGridFile('IgnoreKnm2FSDbinning','ON');
% S.Interp1Grid('Maxm4Sq',38.2^2);
% pHandle{2} =  S.ContourPlot('HoldOn','ON','SavePlot','OFF','Color',rgb('Orange'),'LineStyle','-.');
% HoldOn = 'ON';

% %%
% leg = legend([pHandle{:}],'Total','Stat. only',sprintf('Stat. and {\\itB}_{pt}'));
% PrettyLegendFormat(leg);

%save
% plotname999 = sprintf('%sBudget%.0f_SystBreakdownDiff.png',extractBefore(S.DefPlotName,'Budget'),A.SysBudget);
% print(f999,plotname999,'-dpng','-r300');
% fprintf('save plot to %s \n',plotname999)