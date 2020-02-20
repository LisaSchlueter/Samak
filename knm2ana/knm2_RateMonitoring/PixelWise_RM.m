%% Load KNM2 Period RW2 REAL data for comparison
RunList    = 'KNM1rm';
Knm2AnaArg = {'RunList',RunList,'DataType','Real',...
    'FSDFlag','Sibille0p5eV','ELossFlag','KatrinT2',...
    'AnaFlag','StackPixel','RingMerge','Full','AnaFlag','StackPixel','NonPoissonScaleFactor',1.064};
DataUni_RW2   = MultiRunAnalysis(Knm2AnaArg{:});
DataPSR_RW2   = RingAnalysis('RunAnaObj',DataUni_RW2,'RingList',1:4);

%% Time in days
StartTimeStampDays           = days(DataUni_RW2.SingleRunData.StartTimeStamp-DataUni_RW2.SingleRunData.StartTimeStamp(1));
    
%% Get Activity Correction Per Scan
ScanActCorr   = (DataUni_RW2.SingleRunData.WGTS_MolFrac_TT'+0.5*DataUni_RW2.SingleRunData.WGTS_MolFrac_HT'+0.5*DataUni_RW2.SingleRunData.WGTS_MolFrac_DT')./mean((DataUni_RW2.SingleRunData.WGTS_MolFrac_TT'+0.5*DataUni_RW2.SingleRunData.WGTS_MolFrac_HT'+0.5*DataUni_RW2.SingleRunData.WGTS_MolFrac_DT')).*DataUni_RW2.SingleRunData.WGTS_CD_MolPerCm2'./mean(DataUni_RW2.SingleRunData.WGTS_CD_MolPerCm2);
MeanActivity  = mean((DataUni_RW2.SingleRunData.WGTS_MolFrac_TT'+0.5*DataUni_RW2.SingleRunData.WGTS_MolFrac_HT'+0.5*DataUni_RW2.SingleRunData.WGTS_MolFrac_DT')).*mean(DataUni_RW2.SingleRunData.WGTS_CD_MolPerCm2);

return;

ScanActTable = repmat(ScanActCorr',148,1);
figure(999);set(gcf,'Position',  [100, 100, 1200, 400]);
plot(StartTimeStampDays,ScanActCorr,'LineWidth',3);
xlabel('day');
ylabel('activity corrections');
title(sprintf('Activity Correction - %s - std = %.g',RunList,std(ScanActCorr)));
PrettyFigureFormat    

%% Retreive Pixel-Wise Data - Uniform
DataUni_RW2.ReadSingleRunData;


%% Get 300V sub-scan Time per scan
SStime      = DataUni_RW2.SingleRunData.qUfrac_RM.*DataUni_RW2.SingleRunData.TimeSec;

%% Get qU Correction
qUCorr = (DataUni_RW2.SingleRunData.qU_RM(148,:)); qUCorr=qUCorr-mean(qUCorr);
figure(998);set(gcf,'Position',  [100, 100, 1200, 400]);
plot(StartTimeStampDays,qUCorr,'LineWidth',3);
xlabel('day');
ylabel('qU corrections');
title(sprintf('qU-mean(qU) mV - %s - std = %.g',RunList,std(qUCorr)));
PrettyFigureFormat    

%% Plot Image of Pixel-wise Rate Verus Scan
figure(1)
mapPixelScan = DataUni_RW2.SingleRunData.TBDIS_RM(DataUni_RW2.PixList,:)./ScanActTable(DataUni_RW2.PixList,:)./SStime(DataUni_RW2.PixList,:);
imagesc(mapPixelScan);
%imagesc(DataUni_RW2.SingleRunData.TBDIS_RM(D.PixList,:)./SStime(DataUni_RW2.PixList,:));
colorbar;
xlabel('scan number');
ylabel('golden pixel number');
title(sprintf('%s',RunList));
PrettyFigureFormat
figure(11)
set(gcf,'Position',  [100, 100, 1200, 400]);
stairs(mapPixelScan);
xlabel('golden pixel number');
ylabel('rate per pixel for all scans (cps)');
title(sprintf('%s',RunList));
PrettyFigureFormat
figure(12)
set(gcf,'Position',  [100, 100, 1200, 400]);
stairs(mapPixelScan');
xlabel('scan number');
ylabel('rate per scan for all pixels (cps)');
title(sprintf('%s',RunList));
PrettyFigureFormat
figure(14)
set(gcf,'Position',  [100, 100, 1200, 400]);
stairs(std(mapPixelScan,0,2));
xlabel('golden pixel number');
ylabel('std rate per pixel for all scans (cps)');
title(sprintf('%s',RunList));
PrettyFigureFormat
figure(15)
set(gcf,'Position',  [100, 100, 1200, 400]);
stairs(std(mapPixelScan,0,1));
xlabel('scan number');
ylabel('std rate per scan for all pixels (cps)');
title(sprintf('%s',RunList));
PrettyFigureFormat

%% Retreive Pixel-Wise Data - Ring-wise
for i=1:1:numel(DataPSR_RW2.RingList)
    DataPSR_RW2.MultiObj(i).ReadSingleRunData
end
%% Plot Image of Pixel-wise Rate Verus Scan
figure(2)
for i=1:1:numel(DataPSR_RW2.RingList)
PixelScanMapPSR{i} = DataPSR_RW2.MultiObj(i).SingleRunData.TBDIS_RM(DataPSR_RW2.MultiObj(i).PixList,:)./ScanActTable(DataPSR_RW2.MultiObj(i).PixList,:)./SStime(DataPSR_RW2.MultiObj(i).PixList,:);
subplot(floor(numel(DataPSR_RW2.RingList)/3),floor(numel(DataPSR_RW2.RingList)/3),i)
imagesc((PixelScanMapPSR{i}));
colorbar;
xlabel('scan number');
ylabel('golden pixel number');
title(sprintf('%s - PSR=%.0f - %.2f +- %.2f cps',RunList,i,mean(PixelScanMapPSR{i}(:)),std(PixelScanMapPSR{i}(:))));
PrettyFigureFormat

MeanRatePerPixelPSR(i) = mean(PixelScanMapPSR{i}(:));
StdRatePerPixelPSR(i)  = std(PixelScanMapPSR{i}(:));
end

%% Relative Normalisation WRT PSR1
figure(1000)
set(gcf,'Position',  [100, 100, 1200, 800]);
errorbar(MeanRatePerPixelPSR./MeanRatePerPixelPSR(1),StdRatePerPixelPSR./MeanRatePerPixelPSR./sqrt(numel(DataUni_RW2.RunList)),...
's','MarkerSize',5,'MarkerEdgeColor', rgb('DarkGray'),'LineWidth',3);
xlabel('Pseudo  Ring Number');
ylabel('Normalization wrt PSR1');
title(sprintf('Relative Normalization wrt PSR1 - Per Pixel - %s',RunList));
PrettyFigureFormat
xlim([0.5 9.5])
% xlim([0.5 4.5])
% subplot(2,1,2)
% PixelPSR = [28 36 34 19];
% PixelPSRFactor = PixelPSR./PixelPSR(1);
% errorbar(MeanRatePerPixelPSR./MeanRatePerPixelPSR(1).*PixelPSRFactor,StdRatePerPixelPSR./MeanRatePerPixelPSR.*PixelPSRFactor,...
% 's','MarkerSize',5,'MarkerEdgeColor', rgb('DarkGray'),'LineWidth',3);
% xlabel('Pseudo  Ring Number');
% ylabel('Normalization wrt PSR1');
% title(sprintf('Relative Normalization wrt PSR1 - Per Pixel - %s',RunList));
% PrettyFigureFormat
% xlim([0.5 4.5])

%corrcoef(DataUni_RW2.SingleRunData.EffCorr_RM(DataUni_RW2.PixList,2)',sum(mapPixelScan'))

%% 
figure(10000)
set(gcf,'Position',  [100, 100, 1200, 800]);
stairs(DataUni_RW2.SingleRunData.qU_RM(:,1)-DataUni_RW2.SingleRunData.qU_RM(1,1),'LineWidth',3)
xlabel('pixel')
ylabel('qU(pixel) - qU(pixel0)');
PrettyFigureFormat