% Test SterileAnalysis class
% Lisa, May2020
% plot with neutrino mass as
% nuissance parameter free
% nuissance parameter + pull
% fixed parameter
%% settings for runanalysis
DataType = 'Real';
%%
RunAnaArg = {'RunList','KNM1',...
    'fixPar','mNu E0 Norm Bkg',...
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
%S.InterpMode = 'spline'; %'spline' sometimes causes weird artefacts, but looks smoother than 'lin'
Arg = {'SavePlot','png','BestFit','OFF'};
S.PlotPRL1(Arg{:});
