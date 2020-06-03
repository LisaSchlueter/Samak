% uniform fit on knm2 stacked data
% settings
RunList = 'KNM1';
fixPar = 'E0 Bkg Norm'; % free parameter
DataType = 'Real';
FSDFlag = 'Sibille0p5eV';
NonPoissonScaleFactor = 1.064;
ELossFlag = 'KatrinT2';
AnaFlag = 'StackPixel'; % uniform FPD
range = 40; 
chi2 = 'chi2Stat';
RunAnaArg = {'RunList',RunList,...
             'fixPar',fixPar,...
             'DataType',DataType,...
            'FSDFlag',FSDFlag,...
            'ELossFlag',ELossFlag,...
            'NonPoissonScaleFactor',NonPoissonScaleFactor,...
            'AnaFlag',AnaFlag,...
            'chi2',chi2,...
            'RingMerge','None',...
            'AngularTFFlag','OFF'};

% read data and set up model
A = MultiRunAnalysis(RunAnaArg{:});
A.exclDataStart = A.GetexclDataStart(range);
R = RingAnalysis('RunAnaObj',A,'RingList',1:12);

%%
R.FitRings('SaveResult','OFF','RecomputeFlag','ON');

%%
R.PlotFits('SavePlot','ON','Blind','ON','PlotPar',2,'YLim',[-0.25 0.25])