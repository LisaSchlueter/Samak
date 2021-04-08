% investigate effect of pull terms
% ksn2 calculate chi2 grid search
%% settings that might change
chi2 = 'chi2CMShape';
DataType = 'Twin';
nGridSteps = 30;
range = 40;

mNuSqPullFlag =  [26,12,18:-1:15];

% 12 -> 1.94 eV^2 mainz&troitzk
% 15 -> 1 eV^2
% 16 -> 2 eV^2
% 17 -> 3 eV^2
% 18 -> 0.5 eV^2
% 26 -> 1.1 eV^2 KATRIN KNM-1

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
    'FSDFlag','KNM2_0p1eV',...
    'ELossFlag','KatrinT2A20',...
    'AnaFlag','StackPixel',...
    'chi2',chi2,...
    'NonPoissonScaleFactor',NonPoissonScaleFactor,...
    'FSD_Sigma',sqrt(0.0124+0.0025),...
    'TwinBias_FSDSigma',sqrt(0.0124+0.0025),...
    'TwinBias_Q',18573.7,...
    'BKG_PtSlope',3*1e-06,...
    'TwinBias_BKG_PtSlope',3*1e-06,...
    'DopplerEffectFlag','FSD'};

for i=1:numel(mNuSqPullFlag)
    A = MultiRunAnalysis(RunAnaArg{:},'pullFlag',mNuSqPullFlag(i));
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
    S.GridSearch('ExtmNu4Sq','ON','mNu4SqTestGrid',2);
end