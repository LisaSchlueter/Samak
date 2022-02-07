% uniform fit on knm2 stacked data
% settings
RunList = 'KNM1';
fixPar = 'mNu E0 Bkg Norm'; % free parameter
DataType = 'Real';
FSDFlag = 'Sibille0p5eV';
NonPoissonScaleFactor = 1.064;
ELossFlag = 'KatrinT2';
AnaFlag = 'StackPixel'; % uniform FPD
range = 40; 
chi2 = 'chi2Stat';
RingMerge = 'Full';%None';
RunAnaArg = {'RunList',RunList,...
             'fixPar',fixPar,...
             'DataType',DataType,...
            'FSDFlag',FSDFlag,...
            'ELossFlag',ELossFlag,...
            'NonPoissonScaleFactor',NonPoissonScaleFactor,...
            'AnaFlag',AnaFlag,...
            'chi2',chi2,...
            'RingMerge',RingMerge,...
            'AngularTFFlag','OFF'};

% read data and set up model
A = MultiRunAnalysis(RunAnaArg{:});
A.exclDataStart = A.GetexclDataStart(range);
R = RingAnalysis('RunAnaObj',A,'RingList',A.RingList);

%%
R.FitRings('SaveResult','ON','RecomputeFlag','OFF');

%%
close all
%R.PlotFits('SavePlot','ON','Blind','OFF','PlotPar',1,'PlotMode','Abs','YLim',[-8.5 6.5])
%R.PlotFits('SavePlot','ON','Blind','OFF','PlotPar',2,'PlotMode','Abs','YLim',18573+[+0.25,+1.15])
R.PlotFits('SavePlot','ON','Blind','OFF','PlotPar',3,'PlotMode','Abs','YLim',[2.1 3.1]);%,'YLim',[-8.5 6.5])
%R.PlotFits('SavePlot','ON','Blind','OFF','PlotPar',4,'PlotMode','Rel','YLim',[-1.3e-2 +2.2e-2]);%,'YLim',[-8.5 6.5])


