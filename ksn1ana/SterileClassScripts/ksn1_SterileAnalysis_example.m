% example script how to use SterileAnalysis class
% Lisa, May 2020

% work flow:
% 1. set up a (Multi-)RunAnalysis object "R" with general settings (Runlist, FSD,...)
% 2. set up SterileAnalysis object "S" with "R" as input argument
% 3. from here you can load chi2 grids / plot contours, change some basic settings
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
    'range',65};

S = SterileAnalysis(SterileArg{:});


%% do some thing with the class: 
% 1. change some settings if you want: -> no need to recalc. RunAnalysis object
S.RunAnaObj.DataType = 'Real';
S.range = 40;

%
%S.range = 95;
%S.RunAnaObj.SysBudget = 29;

%% 2 load a chi2 map and find the best fit
S.LoadGridFile('CheckSmallerN','ON','CheckLargerN','ON'); % if CheckSmallerN also look for grid with more/less nGridSteps
S.InterpMode = 'lin'; % waring: if contour is closed, spline interp sometimes sensitive to artefacts! Switch to "lin" in this case
S.Interp1Grid('RecomputeFlag','ON');% interpolate chi2 map -> nicer appearance of all plots. some
S.FindBestFit;
S.CompareBestFitNull;

%% 3. contour plot for some confidence levels

S.RunAnaObj.SysBudget = 24;
S.RunAnaObj.ELossFlag = 'KatrinT2';
S.SysEffect='all';
S.LoadGridFile('CheckSmallerN','ON','CheckLargerN','ON'); % if CheckSmallerN also look for grid with more/less nGridSteps
S.InterpMode = 'spline'; % waring: if contour is closed, spline interp sometimes sensitive to artefacts! Switch to "lin" in this case
S.Interp1Grid('RecomputeFlag','ON');% interpolate chi2 map -> nicer appearance of all plots. some
S.ContourPlot('BestFit','OFF','SavePlot','OFF','CL',[95],'HoldOn','OFF');

S.RunAnaObj.SysBudget = 24;
S.SysEffect='all';
S.RunAnaObj.ELossFlag = 'KatrinT2';
S.LoadGridFile('CheckSmallerN','ON','CheckLargerN','ON'); % if CheckSmallerN also look for grid with more/less nGridSteps
S.InterpMode = 'spline'; % waring: if contour is closed, spline interp sometimes sensitive to artefacts! Switch to "lin" in this case
S.Interp1Grid('RecomputeFlag','ON');% interpolate chi2 map -> nicer appearance of all plots. some
S.ContourPlot('BestFit','OFF','SavePlot','OFF','CL',[95],'HoldOn','ON','LineStyle','-.','Color',rgb('Orange'));
%xlim([1e-03 0.5])
%title(sprintf('%s - %.0f%% C.L.',S.GetPlotTitle,CL),'FontWeight','normal','FontSize',get(gca,'FontSize'));
          
%% 4. grid plot with contour and best fit
S.GridPlot('BestFit','OFF','Contour','ON')

%% 5. contour plot in oscillation parameter space (you can switch on/off foreign contours, all on by default
S.ContourPlotOsci();
 
%% comparisonw ith fitrium (works only for 65 eV & 95 range)
if S.range==95 && strcmp(S.RunAnaObj.DataType,'Real')
    S.InterpMode = 'lin';
else
      S.InterpMode = 'spline';
end
%S.PlotFitriumSamak;

