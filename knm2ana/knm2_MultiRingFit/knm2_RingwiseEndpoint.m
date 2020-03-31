ROIstr = '14keV';


RunList   = 'KNM2_RW2';
Knm2AnaArg = {'RunList',RunList,'fixPar','E0 Bkg Norm','DataType','Real',...
    'FSDFlag','BlindingKNM2','ELossFlag','KatrinT2',...
    'AnaFlag','StackPixel','RingMerge','Full','NonPoissonScaleFactor',1,...
    'MosCorrFlag','OFF',...
    'TwinBias_Q',18573.7,...
    'ROIFlag',ROIstr};
DataUni_RW2   = MultiRunAnalysis(Knm2AnaArg{:});
start = 90;               % fit range in eV below endpoint
stop  = 45;               % fit range in eV below endpoint
DataUni_RW2.exclDataStart = find((DataUni_RW2.ModelObj.qU)>=18574-start,1);
DataUni_RW2.exclDataStop = find((DataUni_RW2.ModelObj.qU)>=18574-stop,1);
DataPSR_RW2   = RingAnalysis('RunAnaObj',DataUni_RW2,'RingList',1);
DataPSR_RW2.MultiObj(1).Fit;
E0Ref = DataPSR_RW2.MultiObj(1).ModelObj.Q_i+DataPSR_RW2.MultiObj(1).FitResult.par(2);


for j = 1:3
    RunList   = ['KNM2_RW' num2str(j)];
    Knm2AnaArg = {'RunList',RunList,'fixPar','E0 Bkg Norm','DataType','Real',...
        'FSDFlag','BlindingKNM2','ELossFlag','KatrinT2',...
        'AnaFlag','StackPixel','RingMerge','Full','NonPoissonScaleFactor',1,...
        'MosCorrFlag','OFF',...
        'TwinBias_Q',18573.7,...
        'ROIFlag',ROIstr};
    DataUni_RW2   = MultiRunAnalysis(Knm2AnaArg{:});
    start = 90;               % fit range in eV below endpoint
    stop  = 45;               % fit range in eV below endpoint
    DataUni_RW2.exclDataStart = find((DataUni_RW2.ModelObj.qU)>=18574-start,1);
    DataUni_RW2.exclDataStop = find((DataUni_RW2.ModelObj.qU)>=18574-stop,1);
    DataPSR_RW2   = RingAnalysis('RunAnaObj',DataUni_RW2,'RingList',1:4);
    for i = 1:DataUni_RW2.nRings
        DataPSR_RW2.MultiObj(i).Fit;
        %DataPSR_RW2.MultiObj(i).FitRunList('Recompute','ON');
        %FitResults = DataPSR_RW2.MultiObj(i).SingleRun_FitResults.chi2Stat;
        E0(i,j) = DataPSR_RW2.MultiObj(i).ModelObj.Q_i+DataPSR_RW2.MultiObj(i).FitResult.par(2);
        %E0Err(i,j) = mean(FitResults.E0Err);
        Diffs_E0(i,j) = E0(i,j)-E0Ref;
    end
end

[PlasmaDrifts, Err_PlasmaDrifts, Offsets, Err_Offsets] = knm2_RingwiseFit('ROI',ROIstr);

fig88 = figure('Renderer','painters');
set(fig88,'units','normalized','pos',[10, 10,1400,1000]);
PlotStyle = { 'o','MarkerSize',8,'MarkerFaceColor',rgb('SkyBlue'),'LineWidth',2};
for i = 1:3
    subplot(3,1,i);
    e1 = errorbar(-Diffs_E0(:,i).*1e3,zeros(4,1),PlotStyle{:});
    hold on;
    e2 = errorbar(-Offsets(:,i),Err_Offsets(:,i),PlotStyle{:});
    xlabel(sprintf('Pseudoring'));
    ylabel(sprintf('mV'));
    leg = legend(sprintf('RW%i {\\Delta}{\\itE}_0',i),sprintf('RW%i {\\Delta}V',i));
    legend boxoff;
    axis([0.5 4.5 -inf inf]);
    hold off;
    PrettyFigureFormat;
end