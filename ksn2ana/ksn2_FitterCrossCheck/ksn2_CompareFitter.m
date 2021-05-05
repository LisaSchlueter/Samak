% ksn2 calculate chi2 grid search
%% settings that might change
chi2 = 'chi2Stat';
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


if strcmp(DataType,'Real')
    S.LoadGridArg = {'mNu4SqTestGrid',5,'IgnoreKnm2FSDbinning','OFF','ExtmNu4Sq','ON'};
else
    S.LoadGridArg = {'mNu4SqTestGrid',5,'IgnoreKnm2FSDbinning','ON','ExtmNu4Sq','OFF'};
end
%%
S.NullHypothesis = 'ON';
S.InterpMode = 'spline';
S.nGridSteps = 35;
S.LoadGridFile('mNu4SqTestGrid',5,'IgnoreKnm2FSDbinning','OFF','ExtmNu4Sq','OFF');
S.Interp1Grid;
p2 = S.ContourPlot('BestFit','ON');

S.nGridSteps = 40;
S.LoadGridFile('mNu4SqTestGrid',5,'IgnoreKnm2FSDbinning','OFF','ExtmNu4Sq','OFF');
S.Interp1Grid;
p1 = S.ContourPlot('BestFit','ON','HoldOn','ON','LineStyle',':','Color',rgb('FireBrick'));

S.nGridSteps = 30;
S.LoadGridFile('mNu4SqTestGrid',5,'IgnoreKnm2FSDbinning','OFF','ExtmNu4Sq','ON');
S.Interp1Grid;
p3 = S.ContourPlot('BestFit','ON','HoldOn','ON','LineStyle','-.','Color',rgb('Orange'));

leg = legend([p3,p2,p1],'30x30','35x35','40x40');
leg.Title.String = 'Grid size';
leg.Title.FontWeight = 'normal';
PrettyLegendFormat(leg);

xlim([5e-03,0.5]);
ylim([1 1600]);
%
savedir = [getenv('SamakPath'),'ksn2ana/ksn2_FitterCrossCheck/plots/'];
MakeDir(savedir);
pltname = sprintf('%sksn2_data_mNuSqFree_GridSize_NH.png',savedir);
print(pltname,'-dpng','-r350');
%%
S.NullHypothesis = 'ON';
S.PlotFitriumSamak('PlotTot','OFF','PlotStat','ON',...
                  'SavePlot','ON','PlotKafit','ON','xLim',[5e-03,0.5])

