% KNM2 uniform fit to 1 (random) single twin run
% Lisa, March 2020
RunAnaArg = {'RunNr',56160,...%56274,...  
    'fixPar','E0 Bkg Norm',...     % free Parameter !!
    'DataType','Real',...              % Real, Twin or Fake
    'FSDFlag','BlindingKNM2',...       % final state distribution (theoretical calculation)
    'ELossFlag','KatrinT2A20',...         % energy loss function     ( different parametrizations available)
    'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
    'chi2','chi2Stat',...              % statistics only
    'NonPoissonScaleFactor',1,...
    'TwinBias_Q',18573.7,...
    'ROIFlag','Default'};

% build object of MultiRunAnalysis class
SR = RunAnalysis(RunAnaArg{:});
SR.exclDataStart = SR.GetexclDataStart(40); % define fit range
% %% enhance statistics for test
% TimeSec = SR.RunData.TimeSec;
% SR.ModelObj.TimeSec = 100*TimeSec;
% SR.RunData.TBDIS = SR.RunData.TBDIS.*100; 

%%
SR.Fit;
%SR.PlotFit;