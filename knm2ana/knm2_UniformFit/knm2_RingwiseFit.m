RunList    = 'KNM2_RW1';
Knm2AnaArg = {'RunList',RunList,'DataType','Real',...
    'FSDFlag','Sibille0p5eV','ELossFlag','KatrinT2',...
    'AnaFlag','StackPixel','RingMerge','Full','AnaFlag','StackPixel','NonPoissonScaleFactor',1};
DataUni_RW2   = MultiRunAnalysis(Knm2AnaArg{:});
range = 40;               % fit range in eV below endpoint        
DataUni_RW2.exclDataStart = DataUni_RW2.GetexclDataStart(range); % find correct data, where to cut spectrum
DataPSR_RW2   = RingAnalysis('RunAnaObj',DataUni_RW2,'RingList',1:4);

for i = 1:DataUni_RW2.nRings
    DataPSR_RW2.MultiObj(i).RMCorrection('saveplot','OFF','pixlist',sprintf('ring%i',i));
    DataPSR_RW2.MultiObj(i).PlotFitRunListCorr('Parameterx','time','Parametery','rate300','Fit','ON','saveplot','OFF','pixlist',sprintf('ring%i',i));
    DataPSR_RW2.MultiObj(i).PlotFitRunListCorr('Parameterx','time','Parametery','rate300','Fit','ON','Detrend','ON');
end