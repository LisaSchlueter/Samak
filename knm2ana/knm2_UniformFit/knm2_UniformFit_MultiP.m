% Uniform fit on KNM2 data
% analyze 3 rear wall periods together
% twins with plasma drift
% March 2020, Lisa

E0OffseteV  = [0,0.1,-0.1]';          % per RW-perid (the same for all rings )
DriftPerDay = [6*1e-03, 0, 6*1e-03]'; % per RW-perid (the same for all rings )

[E0,RectWidth,MultiWeights,RectStd] = ConvertPlasmaDriftDay2E0Run('DriftPerDay',DriftPerDay,...
    'E0OffseteV',E0OffseteV,...
    'E0ref',18573.70,...
    'SanityPlot','OFF');

%%
RunAnaArg = {'RunList','KNM2_Prompt',...  % define run number -> see GetRunList
    'fixPar','mNu E0 Bkg Norm',...         % free Parameter !!
    'DataType','Twin',...              % Real, Twin or Fake
    'FSDFlag','BlindingKNM2',...       % final state distribution (theoretical calculation)
    'ELossFlag','KatrinT2',...         % energy loss function     ( different parametrizations available)
    'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
    'chi2','chi2Stat',...              % statistics only
    'NonPoissonScaleFactor',1,...
    'TwinBias_Q',18573.7,...
    'ROIFlag','Default'};

%% build object of MultiRunAnalysis class
A = MultiRunAnalysis(RunAnaArg{:});

%% modify some parameters in your analysis
range = 40;               % fit range in eV below endpoint        
A.exclDataStart = A.GetexclDataStart(range); % find correct data, where to cut spectrum

%% reference fit
A.Fit;
FitResults_ref = A.FitResult;


%% modiy FSDs
if numel(A.TwinBias_Q)>1
    A.ModelObj.LoadFSD('MultiPos',E0OffseteV,...
        'Sigma',RectStd',...
        'MultiWeights',MultiWeights,...
        'SanityPlot','ON',...
        'BinningFactor',2,...
        'Dist','Gauss');
end
%% Fit 
A.Fit;

%% Display result
%A.PlotFit('saveplot','pdf','FitResultsFlag','OFF');

