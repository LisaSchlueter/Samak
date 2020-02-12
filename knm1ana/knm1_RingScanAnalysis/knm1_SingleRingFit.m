% uniform fit on knm2 stacked data
% settings
RunList = 'KNM1';
fixPar = '5 6 7 8 9 10 11 12'; % fixed parameter
DataType = 'Real';
FSDFlag = 'Sibille0p5eV';
NonPoissonScaleFactor = 1.064;
ELossFlag = 'KatrinT2';
AnaFlag = 'StackPixel'; % uniform FPD
exclDataStart = 14; % 
chi2 = 'chi2Stat';
RunAnaArg = {'RunList',RunList,'fixPar',fixPar,'DataType',DataType,...
            'FSDFlag',FSDFlag,'ELossFlag',ELossFlag,...
            'NonPoissonScaleFactor',NonPoissonScaleFactor,'exclDataStart',exclDataStart,...
            'AnaFlag',AnaFlag,'chi2',chi2,...
            'RingMerge','Full'};

% read data and set up model
A = MultiRunAnalysis(RunAnaArg{:});
R = RingAnalysis('RunAnaObj',A,'RingList',1:4);

%%
R.FitRings('SaveResult','OFF','RecomputeFlag','ON');

%%
R.PlotFits('SavePlot','ON','Blind','ON','PlotPar',4)