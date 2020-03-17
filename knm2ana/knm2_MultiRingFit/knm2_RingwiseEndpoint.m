for j = 1:3
    RunList   = ['KNM2_RW' num2str(j)];
    Knm2AnaArg = {'RunList',RunList,'fixPar','E0 Bkg Norm','DataType','Real',...
        'FSDFlag','BlindingKNM2','ELossFlag','KatrinT2',...
        'AnaFlag','StackPixel','RingMerge','Full','NonPoissonScaleFactor',1,...
        'MosCorrFlag','OFF',...
        'TwinBias_Q',18573.7,...
        'ROIFlag','Default'};
    DataUni_RW2   = MultiRunAnalysis(Knm2AnaArg{:});
    start = 90;               % fit range in eV below endpoint
    stop  = 45;               % fit range in eV below endpoint
    DataUni_RW2.exclDataStart = find((DataUni_RW2.ModelObj.qU)>=18574-start,1);
    DataUni_RW2.exclDataStop = find((DataUni_RW2.ModelObj.qU)>=18574-stop,1);
    DataPSR_RW2   = RingAnalysis('RunAnaObj',DataUni_RW2,'RingList',1:4);
    for i = 1:DataUni_RW2.nRings
        DataPSR_RW2.MultiObj(i).Fit;
        E0(i,j) = DataPSR_RW2.MultiObj(i).ModelObj.Q_i+DataPSR_RW2.MultiObj(i).FitResult.par(2);
        Diffs_E0(i,j) = E0(i,j)-E0(1,1);
    end
end

[PlasmaDrifts, Err_PlasmaDrifts, Offsets, Err_Offsets] = knm2_RingwiseFit('ROI','Default');
fig88 = figure('Renderer','painters');
set(fig88,'units','normalized','pos',[0.1, 0.1,1,0.6]);
PlotStyle = { '-o','MarkerSize',8,'MarkerFaceColor',rgb('SkyBlue'),'LineWidth',2};
e1 = errorbar(Diffs_E0.*1e3,zeros(4,3),PlotStyle{:});
hold on;
e2 = errorbar(Offsets,Err_Offsets,PlotStyle{:});
xlabel(sprintf('Pseudoring'));
ylabel(sprintf('mV'));
leg = legend(sprintf('RW1 {\\Delta}{\\itE}_0'),sprintf('RW2 {\\Delta}{\\itE}_0'),sprintf('RW3 {\\Delta}{\\itE}_0'),sprintf('RW1 {\\Delta}V'),sprintf('RW2 {\\Delta}V'),sprintf('RW3 {\\Delta}V'));
hold off;
PrettyFigureFormat;