% KNM2 Figure skating twins
range = 40;
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
A.exclDataStart = A.GetexclDataStart(range);
A.Fit;
FitResult_ref  = A.FitResult;
%%
TimeSec = zeros(3,1);
TimeSec(1) = sum(A.SingleRunData.TimeSec(1:171));
TimeSec(2) = sum(A.SingleRunData.TimeSec(172:268));
TimeSec(3) = sum(A.SingleRunData.TimeSec(269:361));
MultiWeights = TimeSec./sum(TimeSec);
MultiPos = [E0(1),E0(end-120),E0(end)]';
MultiPosRel = MultiPos-wmean(MultiPos,MultiWeights);
A.ModelObj.LoadFSD('MultiPos',MultiPosRel,'MultiWeight',MultiWeights,'SanityPlot','ON');

A.Fit;
FitResult_imp  = A.FitResult;

