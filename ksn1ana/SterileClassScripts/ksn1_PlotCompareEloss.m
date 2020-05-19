% ksn1 - study influence of energy-loss parametrizations
% Lisa, May 2020

% work flow:

%% 1. RunAnalysis object
RunAnaArg = {'RunList','KNM1',...
    'fixPar','E0 Norm Bkg',...
    'DataType','Real',...
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

R = MultiRunAnalysis(RunAnaArg{:});
R.chi2 = 'chi2CMShape';
%% 2. SterileAnalysis class
SterileArg = {'RunAnaObj',R,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
    'nGridSteps',50,...
    'SmartGrid','OFF',...
    'RecomputeFlag','OFF',...
    'SysEffect','all',...
    'RandMC','OFF',...
    'range',95,...
    'SysBudget',29};

S = SterileAnalysis(SterileArg{:});
%% plot
S.range = 95;
CL = 99;

% new (KNM2) energy-loss
S.RunAnaObj.SysBudget= 29;
S.RunAnaObj.ELossFlag = 'KatrinT2';
S.LoadGridFile('CheckSmallerN','ON','CheckLargerN','ON'); % if CheckSmallerN also look for grid with more/less nGridSteps
S.InterpMode = 'lin'; % waring: if contour is closed, spline interp sometimes sensitive to artefacts! Switch to "lin" in this case
S.Interp1Grid('RecomputeFlag','ON');% interpolate chi2 map -> nicer appearance of all plots. some
S.ContourPlot('BestFit','ON','SavePlot','OFF','CL',CL);

% new (KNM2) energy-loss
S.RunAnaObj.SysBudget=299;
S.RunAnaObj.ELossFlag = 'KatrinT2';%%A20';
S.LoadGridFile('CheckSmallerN','ON','CheckLargerN','ON'); % if CheckSmallerN also look for grid with more/less nGridSteps
S.InterpMode = 'lin'; % waring: if contour is closed, spline interp sometimes sensitive to artefacts! Switch to "lin" in this case
S.Interp1Grid('RecomputeFlag','ON');% interpolate chi2 map -> nicer appearance of all plots. some
S.ContourPlot('BestFit','ON','SavePlot','OFF','CL',CL,'Color',rgb('Orange'),'HoldOn','ON');
%%
% old (KNM1) energy-loss
S.RunAnaObj.SysBudget=24;
S.RunAnaObj.ELossFlag = 'KatrinT2';
S.LoadGridFile('CheckSmallerN','ON','CheckLargerN','ON'); % if CheckSmallerN also look for grid with more/less nGridSteps
S.Interp1Grid('RecomputeFlag','ON');% interpolate chi2 map -> nicer appearance of all plots. some
S.ContourPlot('BestFit','OFF','SavePlot','OFF','CL',CL,'HoldOn','ON','Color',rgb('Orange'));

xlim([1e-03 0.5])
ylim([1 1e4])
leg = legend('KNM2 KATRIN energy-loss function (April 2020)','KNM1 KATRIN energy-loss function',...
    'EdgeColor',rgb('Silver'),'Location','southwest');
title(sprintf('%s - %.0f%% C.L.',S.GetPlotTitle,CL),'FontWeight','normal','FontSize',get(gca,'FontSize'));

plotname = sprintf('%s_ElossComparison_%.2gCL.png',S.DefPlotName,CL);
print(gcf,plotname,'-dpng','-r450');