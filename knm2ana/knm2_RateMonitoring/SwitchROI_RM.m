RunList    = 'KNM2_RW2';
Knm2AnaArg = {'RunList',RunList,'DataType','Real',...
    'FSDFlag','Sibille0p5eV','ELossFlag','KatrinT2',...
    'AnaFlag','StackPixel','RingMerge','Full','AnaFlag','StackPixel',...
    'NonPoissonScaleFactor',1.064,'ROIFlag','14keV'};
DataUni_RW2   = MultiRunAnalysis(Knm2AnaArg{:});


%% Retreive Pixel-Wise Data - Uniform
DataUni_RW2.ReadSingleRunData;

% Sub-scan Time
SStime      = DataUni_RW2.SingleRunData.qUfrac_RM.*DataUni_RW2.SingleRunData.TimeSec;

% ROI 14-32 keV
DataUni_RW2.ROIFlag='14keV'; DataUni_RW2.SetROI;
mapPixelScan_ROI14 = DataUni_RW2.SingleRunData.TBDIS_RM(DataUni_RW2.PixList,:)./SStime(DataUni_RW2.PixList,:);

% ROI 22-32 keV
DataUni_RW2.ROIFlag='Default'; DataUni_RW2.SetROI;
mapPixelScan_ROI22= DataUni_RW2.SingleRunData.TBDIS_RM(DataUni_RW2.PixList,:)./SStime(DataUni_RW2.PixList,:);

%% Image the Difference ROI 14-32 keV - ROI 22-32 keV
figure(8000)
imagesc(mapPixelScan_ROI14-mapPixelScan_ROI22);
colorbar
PrettyFigureFormat