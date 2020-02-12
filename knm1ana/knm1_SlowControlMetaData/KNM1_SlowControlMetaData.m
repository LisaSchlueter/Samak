% Read KNM1 Data
[Real , Twin] = Knm1RealTwin_Create('Real','ON','Twin','ON');

% Plot SC Distributions
Real.PlotSCdistributions;

% Plot RhoD TimeLine
Real.PlotSCdata_RhoD;

% Plot RhoD Error TimeLine
Real.PlotSCdata_RhoDError;

%% Plot Isotopologues TimeLines
Real.PlotSCdata_Isotopologue('Molecule','TT');
Real.PlotSCdata_Isotopologue('Molecule','HT');
Real.PlotSCdata_Isotopologue('Molecule','DT');

%% Plot Isotopologues Uncerrtainty TimeLines
Real.PlotSCdata_IsotopologueError('Molecule','TT');
Real.PlotSCdata_IsotopologueError('Molecule','HT');
Real.PlotSCdata_IsotopologueError('Molecule','DT');

%% Tritium Purity
Real.PlotSCdata_TritiumPurity;

%% Tritium Activity
Real.PlotSCdata_TritiumActivity;
