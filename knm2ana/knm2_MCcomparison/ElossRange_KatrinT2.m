

RunAnaArg = {'RunList','KNM2_Prompt',...  % define run number -> see GetRunList
    'fixPar','mNu E0 Bkg Norm',...         % free Parameter !!
    'DataType','Twin',...              % Real, Twin or Fake
    'FSDFlag','BlindingKNM2',...       % final state distribution (theoretical calculation)
    'ELossFlag','KatrinT2',...         % energy loss function     ( different parametrizations available)
    'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
    'chi2','chi2Stat',...              % statistics only
    'NonPoissonScaleFactor',1,...
    'MosCorrFlag','OFF',...
    'TwinBias_Q',18573.7,...
    'ROIFlag','Default'};

%% build object of MultiRunAnalysis class
A = MultiRunAnalysis(RunAnaArg{:});
TBDIS_ref = A.ModelObj.TBDIS;

A.RunData.TBDIS = TBDIS_ref;
A.Fit;
% 
% 
% % Test influence of energy range of energy-loss function
% A = ref_FakeRun_KNM2_RFcomparison;
% A.recomputeRF = RecomputeFlag;
% [~, ElossFunctions] = A.ComputeELossFunction('E',E);