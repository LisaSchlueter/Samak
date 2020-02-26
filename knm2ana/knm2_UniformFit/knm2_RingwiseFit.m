for j = 1:3
    RunList   = ['KNM2_RW' num2str(j)];
    Knm2AnaArg = {'RunList',RunList,'DataType','Real',...
        'FSDFlag','BlindingKNM2','ELossFlag','KatrinT2',...
        'AnaFlag','StackPixel','RingMerge','Full','NonPoissonScaleFactor',1};
    DataUni_RW2   = MultiRunAnalysis(Knm2AnaArg{:});
    range = 40;               % fit range in eV below endpoint        
    DataUni_RW2.exclDataStart = DataUni_RW2.GetexclDataStart(range); % find correct data, where to cut spectrum
    DataPSR_RW2   = RingAnalysis('RunAnaObj',DataUni_RW2,'RingList',1:4);

    for i = 1:DataUni_RW2.nRings
        DataPSR_RW2.MultiObj(i).Fit;
        DataPSR_RW2.MultiObj(i).PlotFit;
        DataPSR_RW2.MultiObj(i).RMCorrection('saveplot','OFF','pixlist',sprintf('ring%i',i));
        DataPSR_RW2.MultiObj(i).PlotFitRunListCorr('Parameterx','time','Parametery','rate300','Fit','ON','saveplot','OFF','pixlist',sprintf('ring%i',i));
        DataPSR_RW2.MultiObj(i).PlotFitRunListCorr('Parameterx','time','Parametery','rate300','Fit','ON','Detrend','ON');
    end
end