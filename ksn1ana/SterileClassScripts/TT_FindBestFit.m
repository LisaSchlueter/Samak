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



%% 40eV
fprintf(2,'40 eV range\n');
S.RunAnaObj.DataType = 'Real';
S.range = 40;
S.LoadGridFile('CheckSmallerN','ON','CheckLargerN','ON'); % if CheckSmallerN also look for grid with more/less nGridSteps
S.InterpMode = 'spline'; % waring: if contour is closed, spline interp sometimes sensitive to artefacts! Switch to "lin" in this case
S.Interp1Grid('RecomputeFlag','ON');% interpolate chi2 map -> nicer appearance of all plots. some
S.FindBestFit;
S.CompareBestFitNull;
fprintf(2,'\n');

% 65eV
fprintf(2,'65 eV range\n');
S.RunAnaObj.DataType = 'Real';
S.range = 65;
S.LoadGridFile('CheckSmallerN','ON','CheckLargerN','ON'); % if CheckSmallerN also look for grid with more/less nGridSteps
S.InterpMode = 'spline'; % waring: if contour is closed, spline interp sometimes sensitive to artefacts! Switch to "lin" in this case
S.Interp1Grid('RecomputeFlag','ON');% interpolate chi2 map -> nicer appearance of all plots. some
S.FindBestFit;
S.CompareBestFitNull;
fprintf(2,'\n');

% 95eV
fprintf(2,'95 eV range\n');
S.RunAnaObj.DataType = 'Real';
S.range = 95;
S.LoadGridFile('CheckSmallerN','ON','CheckLargerN','ON'); % if CheckSmallerN also look for grid with more/less nGridSteps
S.InterpMode = 'lin'; % waring: if contour is closed, spline interp sometimes sensitive to artefacts! Switch to "lin" in this case
S.Interp1Grid('RecomputeFlag','ON');% interpolate chi2 map -> nicer appearance of all plots. some
S.FindBestFit;
S.CompareBestFitNull;
fprintf(2,'\n');
