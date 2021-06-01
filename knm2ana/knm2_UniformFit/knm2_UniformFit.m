% Uniform fit on KNM2 data
% analyze 3 rear wall periods one after the other
% January 2020, Lisa

RunAnaArg = {'RunList','KNM2_RW1',...  % define run number -> see GetRunList
    'fixPar','E0 Bkg Norm',...         % free Parameter !!
    'DataType','Real',...              % Real, Twin or Fake
    'FSDFlag','BlindingKNM2',...       % final state distribution (theoretical calculation)
    'ELossFlag','KatrinT2',...         % energy loss function     ( different parametrizations available)
    'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
    'chi2','chi2Stat',...           % statistics only
    'NonPoissonScaleFactor',1.112,...
    'MosCorrFlag','OFF',...
    'TwinBias_Q',18573.70,...% 18573.7 = default settings, 18574= const ISX, 18575= same Te for all runs
    'ROIFlag','Default',...    
    'DopplerEffectFlag','FSD',...
    'SysBudget',33};%,...
    %'Twin_SameqUFlag','ON'};

%% build object of MultiRunAnalysis class
D = MultiRunAnalysis(RunAnaArg{:});

%% modify some parameters in your analysis
range = 40;               % fit range in eV below endpoint        
D.exclDataStart = D.GetexclDataStart(range); % find correct data, where to cut spectrum

%% Fit -> fit results are in property: A.FitResult
D.Fit('CATS','OFF');

%% Display result
%D.PlotFit('saveplot','OFF','FitResultsFlag','OFF');

