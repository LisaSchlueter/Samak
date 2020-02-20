% apply FSD correction to uniform fit
% idea: 3 periods with shifted rear wall/plasma potential + each period has a linear drift 


RunAnaArg = {'RunList','KNM2_Prompt',...  % define run number -> see GetRunList
    'fixPar','E0 Bkg Norm',...         % free Parameter !!
    'DataType','Real',...              % Real, Twin or Fake
    'FSDFlag','BlindingKNM2',...       % final state distribution (theoretical calculation)
    'ELossFlag','KatrinT2',...         % energy loss function     ( different parametrizations available)
    'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
    'chi2','chi2Stat',...              % statistics only
    'NonPoissonScaleFactor',1};

% build object of MultiRunAnalysis class
A = MultiRunAnalysis(RunAnaArg{:});
 
%% broaden FSDs
MultiWeights = knm2_RWcombi_GetMultiWeights;
MultiPos     = [-0.1,0,0.1];%3*knm2_RWcombi_GetMultiPos_E0fit;
TimeSec = MultiWeights.*A.RunData.TimeSec;
Slope = [6,2,6]'.*1e-03; % slope in eV/day
RectWidth = ConvertPlasmaDrift2Rect(Slope,TimeSec)';

FSDArg = {'MultiPos',MultiPos,'MultiWeights',MultiWeights,...
    'Sigma',RectWidth,'Dist','Rect','SanityPlot','ON','ZoomPlot','OFF'};
A.ModelObj.LoadFSD(FSDArg{:});

%%


