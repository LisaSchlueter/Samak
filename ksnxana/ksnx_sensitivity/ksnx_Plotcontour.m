% KSNX KATRIN final sensitivity
% Lisa, May 2020
% 750 days
% 130 mcps background
% TDR like systematics
%% settings for runanalysis
DataType = 'Fake';
chi2 = 'chi2Stat';
%   
FakeInitFile = @ref_KSNX_KATRIN_Final;
RunAnaArg = {'RunNr',1,...
    'fixPar','E0 Norm Bkg',...
    'DataType','Fake',...
    'FSDFlag','SibilleFull',...
    'ELossFlag','KatrinT2A20',...
    'AnaFlag','StackPixel',...
    'chi2','chi2Stat',...
    'ROIFlag','Default',...
    'SynchrotronFlag','ON',...
    'AngularTFFlag','ON',...
    'ISCSFlag','Edep',...
    'TwinBias_Q',18573.73,...
    'SysBudget',66,...
    'pullFlag',99,...
    'NonPoissonScaleFactor',1,...
    'FakeInitFile',FakeInitFile,...
    'PixList',1:136};

T = RunAnalysis(RunAnaArg{:});
T.chi2 = chi2;
%% settings sterile class
SterileArg = {'RunAnaObj',T,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
    'nGridSteps',50,...
    'SmartGrid','OFF',...
    'RecomputeFlag','OFF',...
    'SysEffect','all',...
    'RandMC','OFF',...
    'range',40};

S = SterileAnalysis(SterileArg{:});
%%
S.InterpMode = 'Mix';
S.LoadGridFile;
S.Interp1Grid('Maxm4Sq',38^2);
S.ContourPlotOsci;
xlim([5e-03 0.5]);

%%