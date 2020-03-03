for j = 1:3
    RunList   = ['KNM2_RW' num2str(j)];
    Knm2AnaArg = {'RunList',RunList,'fixPar','E0 Bkg Norm','DataType','Real',...
        'FSDFlag','BlindingKNM2','ELossFlag','KatrinT2',...
        'AnaFlag','StackPixel','RingMerge','Full','NonPoissonScaleFactor',1};
    DataUni_RW2   = MultiRunAnalysis(Knm2AnaArg{:});
    start = 90;               % fit start in eV below endpoint
    stop  = 45;               % fit stop  in eV below endpoint
    DataUni_RW2.exclDataStart = find((DataUni_RW2.ModelObj.qU)>=18574-start,1);
    DataUni_RW2.exclDataStop = find((DataUni_RW2.ModelObj.qU)>=18574-stop,1);
    DataPSR_RW2   = RingAnalysis('RunAnaObj',DataUni_RW2,'RingList',1:4);
    for i = 1:DataUni_RW2.nRings
        DataPSR_RW2.MultiObj(i).Fit;
        E0(i,j) = DataPSR_RW2.MultiObj(i).ModelObj.Q_i+DataPSR_RW2.MultiObj(i).FitResult.par(2);
    end
end
plot(E0)
xlabel(sprintf('Pseudoring'));
ylabel(sprintf('{\\itE}_0^{fit}'));
leg = legend('RW1','RW2','RW3');
PrettyFigureFormat;