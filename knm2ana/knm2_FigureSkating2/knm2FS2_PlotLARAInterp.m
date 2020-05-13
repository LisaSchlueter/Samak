%KNM2 sanity plot lara subrun value interpolation

range = 40;

E0 = knm2FS_GetE0Twins('SanityPlot','OFF','Mode','FS2');

RunAnaArg = {'RunList','KNM2_Prompt',... % all KNM2 golden runs
    'fixPar','E0 Bkg Norm',...           % free Parameter !!
    'DataType','Real',...              
    'FSDFlag','BlindingKNM2',...       % final state distribution (theoretical calculation)
    'ELossFlag','KatrinT2A20',...         % energy loss function     ( different parametrizations available)
    'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
    'chi2','chi2Stat',...              % statistics only
    'NonPoissonScaleFactor',1.112,...
    'TwinBias_Q',E0,...
    'ROIFlag','Default',...
    'SysBudget',36,...
    'chi2','chi2Stat'};

MC = MultiRunAnalysis(RunAnaArg{:});
MC.exclDataStart = MC.GetexclDataStart(range);

if numel(E0)>0
    FSDArg = {'SanityPlot','OFF','Sigma',std(E0)};
    MC.ModelObj.LoadFSD(FSDArg{:});
end

%% load subrun-wise LARA values
LiveTimeRuns    = hours(MC.SingleRunData.StartTimeStamp-MC.SingleRunData.StartTimeStamp(1));
LiveTimeSubRuns = sort(reshape(repmat(LiveTimeRuns,[MC.ModelObj.nqU,1])+ ...
                   cumsum((MC.RunData.qUfrac.*MC.RunData.TimeSec)./(MC.nRuns*60*60)),MC.ModelObj.nqU*MC.nRuns,1));
                   
MC.ReadSingleRunData('InterpLARA','OFF');
TT_SubRun_before = reshape(MC.SingleRunData.WGTS_MolFrac_TT_SubRun,MC.ModelObj.nqU*MC.nRuns,1);
HT_SubRun_before = reshape(MC.SingleRunData.WGTS_MolFrac_HT_SubRun,MC.ModelObj.nqU*MC.nRuns,1);
DT_SubRun_before = reshape(MC.SingleRunData.WGTS_MolFrac_DT_SubRun,MC.ModelObj.nqU*MC.nRuns,1);

TT_SubRun_before(TT_SubRun_before==0) = NaN; % ignore in plot
HT_SubRun_before(HT_SubRun_before==0) = NaN; % ignore in plot
DT_SubRun_before(DT_SubRun_before==0) = NaN; % ignore in plot

MC.ReadSingleRunData('InterpLARA','ON');
TT_SubRun_after = reshape(MC.SingleRunData.WGTS_MolFrac_TT_SubRun,MC.ModelObj.nqU*MC.nRuns,1);
HT_SubRun_after = reshape(MC.SingleRunData.WGTS_MolFrac_HT_SubRun,MC.ModelObj.nqU*MC.nRuns,1);
DT_SubRun_after = reshape(MC.SingleRunData.WGTS_MolFrac_DT_SubRun,MC.ModelObj.nqU*MC.nRuns,1);

%% plot
figureHandle = figure('Units','normalized','Position',[0.1,0.1,0.7,0.5]);
NormFacTT = mean(TT_SubRun_after);
NormFacHT = mean(HT_SubRun_after);
NormFacDT = mean(DT_SubRun_after);
plot(LiveTimeSubRuns,TT_SubRun_before-NormFacTT,'.','Color',rgb('Silver'));
hold on;
pTT = plot(LiveTimeSubRuns(isnan(TT_SubRun_before)),TT_SubRun_after(isnan(TT_SubRun_before))-NormFacTT,'.','Color',rgb('IndianRed'));
PrettyFigureFormat('FontSize',22);
xlabel('Time (hours)');
ylabel(sprintf('Rel. concentration'))
plot(LiveTimeSubRuns,HT_SubRun_before-NormFacHT,'.','Color',rgb('Silver'));
pHT = plot(LiveTimeSubRuns(isnan(HT_SubRun_before)),HT_SubRun_after(isnan(HT_SubRun_before))-NormFacHT,'.','Color',rgb('DodgerBlue'));
plot(LiveTimeSubRuns,DT_SubRun_before-NormFacDT,'.','Color',rgb('Silver'));
pDT = plot(LiveTimeSubRuns(isnan(DT_SubRun_before)),DT_SubRun_after(isnan(HT_SubRun_before))-NormFacDT,'.','Color',rgb('GoldenRod'));
leg =legend([pTT,pHT,pDT],sprintf('T_2  -  \\langleT_2\\rangle'),...
   sprintf('HT - \\langleHT\\rangle'),sprintf('DT - \\langleDT\\rangle'),'EdgeColor',rgb('Silver'),...
   'Location','northwest');
ylim([-0.012,0.012])
xlim([-20 max(LiveTimeSubRuns+20)]);
export_fig(gcf,[getenv('SamakPath'),'knm2ana/knm2_systematics/plots/knm2_LARAinterp.pdf']);
