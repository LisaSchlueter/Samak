% ksn2 calculate chi2 grid search
% 4 quadrants
%% settings that might change
chi2 = 'chi2CMShape';
DataType = 'Real';
range = 40;
freePar = 'mNu E0 Norm Bkg';
if contains(freePar,'mNu')
    nGridSteps = 40;
    Extsin2T4 = 'OFF';
else
    nGridSteps = 30;
     Extsin2T4 = 'ON';
end
%% configure RunAnalysis object
if strcmp(chi2,'chi2Stat')
    NonPoissonScaleFactor = 1;
elseif  strcmp(chi2,'chi2CMShape')
    NonPoissonScaleFactor = 1.112;
end
RunAnaArg = {'RunList','KNM2_Prompt',...
    'DataType',DataType,...
    'fixPar',freePar,...%free par
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
%% configure Sterile analysis object
SterileArg = {... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
    'nGridSteps',nGridSteps,...
    'SmartGrid','OFF',...
    'RecomputeFlag','OFF',...
    'SysEffect','all',...
    'RandMC','OFF',...
    'range',range};
CommonArg = {'ExtmNu4Sq','ON','mNu4SqTestGrid',5};
%%
A = MultiRunAnalysis(RunAnaArg{:});
S = SterileAnalysis('RunAnaObj',A,SterileArg{:});
S.GridSearch(CommonArg {:},'Extsin2T4',Extsin2T4);

A = MultiRunAnalysis(RunAnaArg{:});
S = SterileAnalysis('RunAnaObj',A,SterileArg{:});
S.GridSearch(CommonArg {:},'Negsin2T4','ON','Extsin2T4','OFF');

A = MultiRunAnalysis(RunAnaArg{:});
S = SterileAnalysis('RunAnaObj',A,SterileArg{:});
S.GridSearch(CommonArg {:},'NegmNu4Sq','ON','Extsin2T4',Extsin2T4);

A = MultiRunAnalysis(RunAnaArg{:});
S = SterileAnalysis('RunAnaObj',A,SterileArg{:});
S.GridSearch(CommonArg {:},'Negsin2T4','ON','NegmNu4Sq','ON','Extsin2T4','OFF');

%% NE quadrant, but extended m4sq and sin2t4
% A = MultiRunAnalysis(RunAnaArg{:});
% S = SterileAnalysis('RunAnaObj',A,SterileArg{:});
% S.GridSearch('ExtmNu4Sq','0.01','mNu4SqTestGrid',2,'Extsin2T4','ON');
% merge quadrant files: tbd
