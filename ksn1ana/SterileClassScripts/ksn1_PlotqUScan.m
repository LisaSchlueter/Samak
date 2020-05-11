% Test SterileAnalysis class
% Lisa, May2020
% plot qU Scan
% fixed parameter
%% settings for runanalysis
DataType = 'Real';
%%
RunAnaArg = {'RunList','KNM1',...
    'fixPar','E0 Norm Bkg',...
    'DataType',DataType,...
    'FSDFlag','SibilleFull',...
    'ELossFlag','KatrinT2',...
    'AnaFlag','StackPixel',...
    'chi2','chi2Stat',...
    'ROIFlag','Default',...
    'SynchrotronFlag','ON',...
    'AngularTFFlag','OFF',...
    'ISCSFlag','Edep',...
    'TwinBias_Q',18573.73,...
    'SysBudget',24,...
    'pullFlag',99,...
    'NonPoissonScaleFactor',1};

T = MultiRunAnalysis(RunAnaArg{:});
T.chi2 = 'chi2CMShape';
%% settings sterile class
SterileArg = {'RunAnaObj',T,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
    'nGridSteps',50,...
    'SmartGrid','OFF',...
    'RecomputeFlag','OFF',...
    'SysEffect','all',...
    'RandMC','OFF',...
    'range',65};

S = SterileAnalysis(SterileArg{:});
%%
%%
switch S.RunAnaObj.DataType
    case 'Real'
        S.InterpMode = 'lin';
        Arg = {'SavePlot','png','Ranges',[95:-5:45,41,40]};
    case 'Twin'
        S.InterpMode = 'spline';
        Arg = {'SavePlot','png','Ranges',[95:-5:45,41,40]};%,41,40]};
end
S.PlotqUScan(Arg{:});

