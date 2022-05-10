% ksn2 90 eV range with free mNu
% remove not converged grid point 
% and plot contour
%% settings that might change
freePar = 'mNu E0 Norm Bkg';
chi2 = 'chi2CMShape';
DataType = 'Real';
nGridSteps = 45;
range = 90;

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
A = MultiRunAnalysis(RunAnaArg{:});
%% configure Sterile analysis object
SterileArg = {... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
    'nGridSteps',nGridSteps,...
    'SmartGrid','OFF',...
    'RecomputeFlag','OFF',...
    'SysEffect','all',...
    'RandMC','OFF',...
    'range',range,...
    'LoadGridArg',{'ExtmNu4Sq','OFF','mNu4SqTestGrid',5},...
    'InterpMode','lin'};
S = SterileAnalysis('RunAnaObj',A,SterileArg{:});

%%
S.InterpMode = 'spline';
S.LoadGridFile(S.LoadGridArg{:});
S.Interp1Grid;
S.GridPlot('Contour','ON','BestFit','ON');
%%
S.InterpMode = 'lin';
S.Interp1Grid('nInter',70)
S.GridPlot('Contour','ON');
chi2_ref = S.chi2_ref;

S.sin2T4 = S.sin2T4';
S.mNu4Sq = S.mNu4Sq';
%S.chi2 = S.chi2';
S.InterpMode = 'spline';
S.Interp1Grid('nInter',1000)
S.chi2_ref = chi2_ref;
S.GridPlot('Contour','ON');

%%

%  
% S.mNu4Sq(:,13) = [];
% S.sin2T4(:,13) = [];
% S.chi2(13,:) = [];

S.Interp1Grid;
S.chi2_ref = 27;
S.ConfLevel = 99;
S.GridPlot('BestFit','ON');





