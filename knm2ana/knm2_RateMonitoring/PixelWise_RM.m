%% This macro get the rate per pixel per scan
%% For a given runlist
%% Avergage rate / period per pseudo-ring
%% Compute Rate Variaton Per Pseud-ring
%% Convert into ROI
%%
%% Allow for KNM1 Calibration
%% Allow to Change ROI
%%
%RunList    = 'KNM2_RW2';
RunList    = 'KNM1rm';
KNM1CorFlag    = 'OFF';
HVdriftCorFlag = 'OFF';
SlopeCPSMeV = 6.3032;

savedir = [getenv('SamakPath'),'knm2ana/knm2_RateMonitoring/results/'];
savename = sprintf('PixelWise_RM_%s_KNM1CorFlag%s_HVdriftCorFlag%s.mat',RunList,KNM1CorFlag,HVdriftCorFlag);


switch RunList
    case 'KNM1rm'
        SlopeCPSMeV = 0.85;
end

    
%% KNM1 with Calibration - Divide Rates by
switch KNM1CorFlag
    case 'OFF'
        KNM1correction  = [ 1.0000    1          1         1     ];
    case 'ON'
        KNM1correction  = [ 1.0000    0.9992    0.9975    0.9961];
end


%% Load KNM2 Period RWX
Knm2AnaArg = {'RunList',RunList,'DataType','Real',...
    'FSDFlag','BlindingKNM2','ELossFlag','KatrinT2',...
    'AnaFlag','StackPixel','RingMerge','Full','AnaFlag','StackPixel',...
    'NonPoissonScaleFactor',1.064,'ROIFlag','14keV'};
DataUni_RW2   = MultiRunAnalysis(Knm2AnaArg{:});
DataPSR_RW2   = RingAnalysis('RunAnaObj',DataUni_RW2,'RingList',1:4);

%% Time in days
StartTimeStampDays           = days(DataUni_RW2.SingleRunData.StartTimeStamp-DataUni_RW2.SingleRunData.StartTimeStamp(1));

%% HV Drift Correction
FirstDayPeriod1 = datetime('02-Oct-2019 14:52:19');
TimeLineDaysFirstDayPeriod1 = days(DataUni_RW2.SingleRunData.StartTimeStamp-FirstDayPeriod1);
switch HVdriftCorFlag
    case 'OFF'
        HVdrift = repmat((0*(TimeLineDaysFirstDayPeriod1) * SlopeCPSMeV*1e-3),numel(DataUni_RW2.PixList),1);
    case 'ON'
        HVdrift = repmat((1.52*(TimeLineDaysFirstDayPeriod1) * SlopeCPSMeV*1e-3),numel(DataUni_RW2.PixList),1);
end
figure(7000)
stairs(DataUni_RW2.SingleRunData.StartTimeStamp,HVdrift(1,:),'LineWidth',3);
xlabel('Date');
ylabel('HV-drift \Delta Rate per Pixel (cps)');
title(sprintf('HV-drift - %s',RunList));
PrettyFigureFormat


%% Get Activity Correction Per Scan Per Period RWX
ScanActCorr   = (DataUni_RW2.SingleRunData.WGTS_MolFrac_TT'+0.5*DataUni_RW2.SingleRunData.WGTS_MolFrac_HT'+0.5*DataUni_RW2.SingleRunData.WGTS_MolFrac_DT')./mean((DataUni_RW2.SingleRunData.WGTS_MolFrac_TT'+0.5*DataUni_RW2.SingleRunData.WGTS_MolFrac_HT'+0.5*DataUni_RW2.SingleRunData.WGTS_MolFrac_DT')).*DataUni_RW2.SingleRunData.WGTS_CD_MolPerCm2'./mean(DataUni_RW2.SingleRunData.WGTS_CD_MolPerCm2);
MeanActivity  = mean((DataUni_RW2.SingleRunData.WGTS_MolFrac_TT'+0.5*DataUni_RW2.SingleRunData.WGTS_MolFrac_HT'+0.5*DataUni_RW2.SingleRunData.WGTS_MolFrac_DT')).*mean(DataUni_RW2.SingleRunData.WGTS_CD_MolPerCm2);

ScanActTable = repmat(ScanActCorr',148,1);
figure(999);set(gcf,'Position',  [100, 100, 1200, 400]);
plot(StartTimeStampDays,ScanActCorr,'LineWidth',3);
xlabel('day');
ylabel('activity corrections');
title(sprintf('Activity Correction - %s - std = %.g',RunList,std(ScanActCorr)));
PrettyFigureFormat

%% Retreive Pixel-Wise Data - Uniform - Change ROI
DataUni_RW2.ReadSingleRunData;
DataUni_RW2.ROIFlag='14keV';
DataUni_RW2.SetROI;

%% Get 300V sub-scan Time per scan
SStime      = DataUni_RW2.SingleRunData.qUfrac_RM.*DataUni_RW2.SingleRunData.TimeSec;

%% Get qU Correction
qUCorr    = (DataUni_RW2.SingleRunData.qU_RM(148,:));
qUmean    = mean(qUCorr);
qUCorr    = qUCorr-mean(qUCorr);
qUCorrCPS = -qUCorr .* SlopeCPSMeV;
qUCorrMap = repmat(qUCorrCPS,numel(DataUni_RW2.PixList),1);

figure(998);set(gcf,'Position',  [100, 100, 1200, 400]);
subplot(2,1,1)
plot(StartTimeStampDays,qUCorr,'LineWidth',3);
xlabel('day');
ylabel('qU corrections');
title(sprintf('qU-mean(qU) mV - %s - std = %.g',RunList,std(qUCorr)));
PrettyFigureFormat
subplot(2,1,2)
plot(StartTimeStampDays,qUCorrCPS,'LineWidth',3);
xlabel('day');
ylabel('cps per pixel');
title(sprintf('qU-mean(qU) mV - %s - std = %.g',RunList,std(qUCorr)));
PrettyFigureFormat

%% Plot Image of Pixel-wise Rate Verus Scan
mapPixelScan = (DataUni_RW2.SingleRunData.TBDIS_RM(DataUni_RW2.PixList,:)...
    ./SStime(DataUni_RW2.PixList,:)+HVdrift+qUCorrMap)...
    ./ScanActTable(DataUni_RW2.PixList,:);

figure(1)
imagesc(mapPixelScan);
colorbar;
xlabel('scan number');
ylabel('golden pixel number');
title(sprintf('%s',RunList));
PrettyFigureFormat

DEBUGPlot = 'OFF';
switch DEBUGPlot
    case 'ON'
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
end

%% Retreive Pixel-Wise Data - Ring-wise
for i=1:1:numel(DataPSR_RW2.RingList)
    DataPSR_RW2.MultiObj(i).ReadSingleRunData
    DataPSR_RW2.MultiObj(i).ROIFlag='14keV';
    DataPSR_RW2.MultiObj(i).SetROI;
    % HV Drift Correction
    switch HVdriftCorFlag
        case 'OFF'
            HVdriftPSR{i} = repmat((0*(TimeLineDaysFirstDayPeriod1) * SlopeCPSMeV*1e-3),numel(DataPSR_RW2.MultiObj(i).PixList),1);
        case 'ON'
            HVdriftPSR{i} = repmat((1.52*(TimeLineDaysFirstDayPeriod1) * SlopeCPSMeV*1e-3),numel(DataPSR_RW2.MultiObj(i).PixList),1);
    end
    % Regular qU Correction
    qUCorrMapPSR{i} = repmat(qUCorrCPS,numel(DataPSR_RW2.MultiObj(i).PixList),1);

end
%% Plot Image of Pixel-wise Rate Verus Scan
figure(2)
for i=1:1:numel(DataPSR_RW2.RingList)
    
PixelScanMapPSR{i} = (HVdriftPSR{i} + ...
    (qUCorrMapPSR{i} + ...
    (DataPSR_RW2.MultiObj(i).SingleRunData.TBDIS_RM(DataPSR_RW2.MultiObj(i).PixList,:))...
    ./SStime(DataPSR_RW2.MultiObj(i).PixList,:))./ScanActTable(DataPSR_RW2.MultiObj(i).PixList,:))...
    ./KNM1correction(i);

subplot(floor(numel(DataPSR_RW2.RingList)/2),floor(numel(DataPSR_RW2.RingList)/2),i)
imagesc((PixelScanMapPSR{i}));
colorbar;
xlabel('scan number');
ylabel('golden pixel number');
title(sprintf('%s - PSR=%.0f - %.2f +- %.2f cps',...
    RunList,i,...
    mean(PixelScanMapPSR{i}(:)),std(PixelScanMapPSR{i}(:))./sqrt(numel(DataPSR_RW2.MultiObj(i).PixList))));
PrettyFigureFormat

MeanRatePerPixelPSR(i) = mean(PixelScanMapPSR{i}(:));
StdRatePerPixelPSR(i)  = std(PixelScanMapPSR{i}(:));
end

%% Relative Normalisation / mV WRT PSR1
figure(1000)
set(gcf,'Position',  [100, 100, 1200, 800]);
subplot(2,1,1)
RateWrtPSR1cps      = MeanRatePerPixelPSR./MeanRatePerPixelPSR(1);
RateErrorWrtPSR1cps = StdRatePerPixelPSR./MeanRatePerPixelPSR./sqrt(numel(DataUni_RW2.RunList));
errorbar(RateWrtPSR1cps,RateErrorWrtPSR1cps,'s','MarkerSize',15,'MarkerEdgeColor', rgb('IndianRed'),'LineWidth',3);
xlabel('Pseudo  Ring Number');
ylabel('Normalization wrt PSR1');
title(sprintf('Relative Normalization wrt PSR1 - Per Pixel - %s',RunList));
PrettyFigureFormat
xlim([0.5 4.5]); xticks([1 2 3 4]);
subplot(2,1,2)
PotentialWrtPSR1mV       = -(MeanRatePerPixelPSR-MeanRatePerPixelPSR(1))./SlopeCPSMeV*1e3;
PotentialErrorWrtPSR1mV  = StdRatePerPixelPSR./sqrt(numel(DataUni_RW2.RunList))./KNM1correction(i)./SlopeCPSMeV*1e3;
errorbar(PotentialWrtPSR1mV,PotentialErrorWrtPSR1mV,'s','MarkerSize',15,'MarkerEdgeColor', rgb('IndianRed'),'LineWidth',3);
xlabel('Pseudo  Ring Number');
ylabel('mV PP shift wrt PSR1');
title(sprintf('Relative Plasma Potential Bound wrt PSR1 - Per Pixel - %s',RunList));
PrettyFigureFormat
xlim([0.5 4.5]);xticks([1 2 3 4]);

table(RateWrtPSR1cps)
table(RateErrorWrtPSR1cps)
table(PotentialWrtPSR1mV)
table(PotentialErrorWrtPSR1mV)
table(MeanRatePerPixelPSR)
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
% figure(10000)
% set(gcf,'Position',  [100, 100, 1200, 800]);
% stairs(DataUni_RW2.SingleRunData.qU_RM(:,1)-DataUni_RW2.SingleRunData.qU_RM(1,1),'LineWidth',3)
% xlabel('pixel')
% ylabel('qU(pixel) - qU(pixel0)');
% PrettyFigureFormat
qU_RM = DataUni_RW2.RunData.qU_RM;
PixList =  DataUni_RW2.PixList;
PixelMap = cell2mat(PixelScanMapPSR');
RingPixList = DataUni_RW2.RingPixList;
save([savedir,savename],'PixelScanMapPSR','qUCorrMapPSR',...
    'PixelMap','qU_RM','PixList','RingPixList');
