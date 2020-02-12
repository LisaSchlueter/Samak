% Read KNM1 Data
[Real , Twin] = Knm1RealTwin_Create('Real','ON','Twin','ON');


% Define the Range at -40 eVs
Real.exclDataStart = 14;
% Retreive Fit Results
Real.FitRunList('Recompute','ON');

% Plot Run-Wise Fit Distributions
Real.PlotFitResultsDistributions
Real.PlotFitRunList

% Stability Plots
FitChi2minFigures;
FitEndpointStability;
FitBackgroundStability;

% Define the Range at -94 eVs
Real.exclDataStart = 2;
% Retreive Fit Results
Real.FitRunList('Recompute','ON');

% Plot Run-Wise Fit Distributions
Real.PlotFitResultsDistributions
Real.PlotFitRunList
