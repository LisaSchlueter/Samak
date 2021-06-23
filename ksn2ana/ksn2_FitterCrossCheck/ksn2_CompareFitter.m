% ksn2 calculate chi2 grid search
%% settings that might change
chi2 = 'chi2Stat';%CMShape';
DataType = 'Twin';
nGridSteps = 50;
range = 40;

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
    'range',range};
S = SterileAnalysis(SterileArg{:});

%%
if strcmp(DataType,'Real')
    S.LoadGridArg = {'mNu4SqTestGrid',5,'IgnoreKnm2FSDbinning','OFF','ExtmNu4Sq','OFF'};
else
    S.LoadGridArg = {'mNu4SqTestGrid',5,'IgnoreKnm2FSDbinning','ON','ExtmNu4Sq','OFF'};
end


%%
 S.LoadGridArg = {'mNu4SqTestGrid',5,'IgnoreKnm2FSDbinning','ON','ExtmNu4Sq','OFF'};
S.nGridSteps = 40;
S.NullHypothesis = 'OFF';
S.PlotFitriumSamak('PlotTot','OFF','PlotStat','ON',...
                  'SavePlot','png','PlotKafit','OFF','xLim',[2e-03,0.5])

          
S.nGridSteps = 30;
S.LoadGridFile(S.LoadGridArg{:});
S.Interp1Grid;
S.ContourPlot('HoldOn','ON','LineStyle',':');