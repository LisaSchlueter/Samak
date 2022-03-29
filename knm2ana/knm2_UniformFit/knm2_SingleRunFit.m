% neutrino-mass sensitivity for single scan
RunAnaArg = {'RunNr',57128,...%,...  % define run number -> see GetRunList
    'fixPar','mNu E0 Bkg Norm',...         % free Parameter !!
    'DataType','Twin',...              % Real, Twin or Fake
    'FSDFlag','KNM2_0p1eV',...       % final state distribution (theoretical calculation)
    'ELossFlag','KatrinT2A20',...         % energy loss function     ( different parametrizations available)
    'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
    'chi2','chi2Stat',...              % statistics only
    'NonPoissonScaleFactor',1.112,...
    'BKG_PtSlope',3*1e-06,...
    'TwinBias_BKG_PtSlope',3*1e-06,...
    'TwinBias_Q',18573.7,...
    'FSD_Sigma',sqrt(0.0124+0.0025),...
    'TwinBias_FSDSigma',sqrt(0.0124+0.0025),...
    'DopplerEffectFlag','FSD'};

%run 56278:  sigma(m^2) = 6.5 eV^2 (twin)
%% build object of MultiRunAnalysis class
A = RunAnalysis(RunAnaArg{:});

%% modify some parameters in your analysis
range = 40; % fit range in eV below endpoint
A.exclDataStart = A.GetexclDataStart(range); % find correct data, where to cut spectrum

%% Fit -> fit results are in property: A.FitResult
A.Fit;

%% Display result
%A.PlotFit('saveplot','OFF');
