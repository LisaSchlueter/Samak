% ksn2 calculate chi2 grid search
%% settings that might change
chi2 = 'chi2CMShape';
DataType = 'Real';
nGridSteps = 30;
range = 40;

%% configure RunAnalysis object
if strcmp(chi2,'chi2Stat')
    NonPoissonScaleFactor = 1;
elseif  strcmp(chi2,'chi2CMShape')
    NonPoissonScaleFactor = 1.112;
end

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
    'LoadGridArg',{'mNu4SqTestGrid',5,'ExtmNu4Sq','ON'}};

%%
S = SterileAnalysis(SterileArg{:});
%S.GridSearch('ExtmNu4Sq','ON');
%%
A.chi2 = 'chi2CMShape';
S.SetNPfactor;
A.fixPar ='fix 5 ;fix 6 ;fix 7 ;fix 8 ;fix 9 ;fix 10 ;fix 11 ;fix 12 ;fix 13 ;fix 14 ;fix 15 ;fix 16 ;fix 17 ;';
S.LoadGridFile('ExtmNu4Sq','OFF');

if strcmp(A.DataType,'Real')
    S.InterpMode = 'spline';
    BF = 'ON';
else
    S.InterpMode = 'spline';
   BF = 'OFF';
end
S.Interp1Grid('RecomputeFlag','ON','Maxm4Sq',39^2);
S.GridPlot('Contour','ON','BestFit',BF,'SavePlot','OFF','CL',95)
%S.ContourPlot('BestFit','OFF','CL',95)
% S.PlotStatandSys('SavePlot','png')
%S.PlotmNuSqOverview('PullmNuSq','OFF','SavePlot','png')