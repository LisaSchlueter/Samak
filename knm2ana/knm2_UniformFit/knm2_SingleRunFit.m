RunAnaArg = {'RunNr',56278,...  % define run number -> see GetRunList
    'fixPar','E0 Bkg Norm',...         % free Parameter !!
    'DataType','Real',...              % Real, Twin or Fake
    'FSDFlag','BlindingKNM2',...       % final state distribution (theoretical calculation)
    'ELossFlag','KatrinT2',...         % energy loss function     ( different parametrizations available)
    'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
    'chi2','chi2Stat',...              % statistics only
    'NonPoissonScaleFactor',1,...
    'MosCorrFlag','ON'};        

%% build object of MultiRunAnalysis class
A = RunAnalysis(RunAnaArg{:});

%% modify some parameters in your analysis
range = 40; % fit range in eV below endpoint
A.exclDataStart = A.GetexclDataStart(range); % find correct data, where to cut spectrum

%% Fit -> fit results are in property: A.FitResult
A.Fit;

%% Display result
A.PlotFit('saveplot','OFF');
