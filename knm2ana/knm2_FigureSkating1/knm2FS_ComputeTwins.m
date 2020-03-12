% KNM2 Figure skating twins

E0 = knm2FS_GetE0Twins('SanityPlot','OFF');

RunAnaArg = {'RunList','KNM2_Prompt',... % all KNM2 golden runs
    'fixPar','mNu E0 Bkg Norm',...           % free Parameter !!
    'DataType','Twin',...              
    'FSDFlag','BlindingKNM2',...       % final state distribution (theoretical calculation)
    'ELossFlag','KatrinT2',...         % energy loss function     ( different parametrizations available)
    'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
    'chi2','chi2Stat',...              % statistics only
    'NonPoissonScaleFactor',1,...
    'TwinBias_Q',E0,...
    'ROIFlag','14keV'}; 

%% build object of MultiRunAnalysis class
A = MultiRunAnalysis(RunAnaArg{:});
%%
A.exclDataStart = A.GetexclDataStart(40);
A.Fit;


