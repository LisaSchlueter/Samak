% ksn2
% m2 nuisance parameter unconstrained
% closer look at area [0.5,1] and compare
%% settings that might change
chi2 = 'chi2Stat';
DataType = 'Twin';
nGridSteps = 30;
range = 40;
Mode = 'Compute';
%% configure RunAnalysis object
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

S = SterileAnalysis(SterileArg{:});

if strcmp(A.DataSet,'Twin')
    ExtmNu4Sq = 'OFF';
else
    ExtmNu4Sq = 'ON';
end
%%
switch Mode
    case 'Compute'
        S.GridSearch('ExtmNu4Sq',ExtmNu4Sq,'mNu4SqTestGrid',5.5); 
        A = MultiRunAnalysis(RunAnaArg{:});
        S = SterileAnalysis(SterileArg{:});    
        S.GridSearch('ExtmNu4Sq',ExtmNu4Sq,'mNu4SqTestGrid',5.5);
    case 'Display'
end
