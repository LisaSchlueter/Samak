RunAnaArg = {'RunList','KNM2_Prompt',...  % define run number -> see GetRunList
    'fixPar','E0 Bkg Norm',...         % free Parameter !!
    'DataType','Real',...              % Real, Twin or Fake
    'FSDFlag','BlindingKNM2',...       % final state distribution (theoretical calculation)
    'ELossFlag','KatrinT2',...         % energy loss function     ( different parametrizations available)
    'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
    'chi2','chi2Stat',...              % statistics only
    'NonPoissonScaleFactor',1,...
    'MosCorrFlag','OFF',...
    'TwinBias_Q',18573.73,...% 18573.7 = default settings, 18574= const ISX, 18575= same Te for all runs
    'ROIFlag','14keV',...    % 18573.73 == new RF binning without interp
    'DopplerEffectFlag','FSD'};%,...


% build object of MultiRunAnalysis class
D = MultiRunAnalysis(RunAnaArg{:});
D.ReadSingleRunData;
BkgRatePixel = NaN.*zeros(148,1);
BkgRatePixel(D.PixList) = mean(squeeze(sum(D.SingleRunData.TBDIS(D.RunData.qU>18574,:,D.PixList),2))./...
    (D.RunData.qUfrac(D.RunData.qU>18574).*D.RunData.TimeSec));
%% relative background rate
[plotHandle, cbHandle] = FPDViewer(100.*(BkgRatePixel-nanmean(BkgRatePixel))./nanmean(BkgRatePixel));
cbHandle.Label.String = 'Rel. background rate (%)';
cbHandle.Label.FontSize = get(gca,'FontSize')+4;
savedir= [getenv('SamakPath'),'knm2ana/knm2_background/plots/'];
MakeDir(savedir);
savename = sprintf('%sknm2_FPDrelBackground.png',savedir);
print(plotHandle,savename,'-dpng','-r400');


%% absolute background rate
[plotHandle, cbHandle] = FPDViewer(1e3.*BkgRatePixel);
cbHandle.Label.String = 'Background rate (mcps)';
cbHandle.Label.FontSize = get(gca,'FontSize')+4;
savedir= [getenv('SamakPath'),'knm2ana/knm2_background/plots/'];
MakeDir(savedir);
savename = sprintf('%sknm2_FPDAbsBackground.png',savedir);
print(plotHandle,savename,'-dpng','-r400');