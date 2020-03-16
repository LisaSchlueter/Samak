% compute and fit twins with a constant ISX

RunAnaArg = {'RunList','KNM2_Prompt',...  % define run number -> see GetRunList
    'fixPar','mNu E0 Bkg Norm',...         % free Parameter !!
    'DataType','Twin',...              % Real, Twin or Fake
    'FSDFlag','BlindingKNM2',...       % final state distribution (theoretical calculation)
    'ELossFlag','KatrinT2',...         % energy loss function     ( different parametrizations available)
    'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
    'chi2','chi2Stat',...              % statistics only
    'NonPoissonScaleFactor',1,...
    'MosCorrFlag','OFF',...
    'TwinBias_Q',18574,...% different E0 to keep track that this is with const. ISX
    'ROIFlag','14keV',...
    'DopplerEffectFlag','FSD',...
    'ISCS','Theory'};
   

%% build object of MultiRunAnalysis class
D = MultiRunAnalysis(RunAnaArg{:});

%% modify some parameters in your analysis
range = 40;               % fit range in eV below endpoint        
D.exclDataStart = D.GetexclDataStart(range); % find correct data, where to cut spectrum

%% Fit -> fit results are in property: A.FitResult
D.Fit;

%% Display result
D.PlotFit('saveplot','OFF','FitResultsFlag','OFF');

