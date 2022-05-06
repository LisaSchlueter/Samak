% ksn2 calculate chi2 grid search
%% settings that might change
DataType              = 'Real';
range                 = 40;
chi2                  = 'chi2CMShape';
freePar               = 'mNu E0 Norm Bkg';

if contains(freePar,'mNu')
    nGridSteps = 40;
    Extsin2T4 = 'OFF';
else
    nGridSteps = 30;
     Extsin2T4 = 'ON';
end

% configure RunAnalysis object
if strcmp(chi2,'chi2Stat')
    NonPoissonScaleFactor = 1;
elseif  strcmp(chi2,'chi2CMShape')
    NonPoissonScaleFactor = 1.064;
end
RunAnaArg = {'RunList','KNM1',...
    'chi2',chi2,...
    'DataType',DataType,...
    'fixPar',freePar,...
    'NonPoissonScaleFactor',NonPoissonScaleFactor,...
    'SysBudget',200,...
    'minuitOpt','min ; minos',...
    'FSDFlag','KNM2_0p1eV',...
    'ELossFlag','KatrinT2A20',...
    'AngularTFFlag','ON',...
    'SynchrotronFlag','ON',...
    'RadiativeFlag','ON',...
    'DopplerEffectFlag','FSD',...
    'BKG_PtSlope',-2.2*1e-06};

%% configure Sterile analysis object
SterileArg = {... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
    'nGridSteps',nGridSteps,...
    'SmartGrid','OFF',...
    'RecomputeFlag','OFF',...
    'SysEffect','all',...
    'RandMC','OFF',...
    'range',range};

CommonArg = {'ExtmNu4Sq','ON','mNu4SqTestGrid',2};
%
A = MultiRunAnalysis(RunAnaArg{:});
S = SterileAnalysis('RunAnaObj',A,SterileArg{:});

%%
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


%%
% S.LoadGridArg = CommonArg;
% S.PlotQuadrant('SavePlot','OFF');