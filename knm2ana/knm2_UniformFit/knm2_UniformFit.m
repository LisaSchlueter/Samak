% Uniform fit on KNM2 data
% analyze 3 rear wall periods one after the other
% January 2020, Lisa

RunAnaArg = {'RunList','KNM2_Prompt',...  % define run number -> see GetRunList
    'fixPar','mNu E0 Bkg Norm',...         % free Parameter !!
    'DataType','Twin',...              % Real, Twin or Fake
    'FSDFlag','BlindingKNM2',...       % final state distribution (theoretical calculation)
    'ELossFlag','KatrinT2',...         % energy loss function     ( different parametrizations available)
    'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
    'chi2','chi2Stat',...              % statistics only
    'NonPoissonScaleFactor',1,...
    'MosCorrFlag','OFF',...
    'TwinBias_Q',18575,...
    'ROIFlag','14keV',...
    'DopplerEffectFlag','FSD'};%,...
    %'Twin_SameqUFlag','ON'};

%% build object of MultiRunAnalysis class
D = MultiRunAnalysis(RunAnaArg{:});

%% modify some parameters in your analysis
range = 40;               % fit range in eV below endpoint        
D.exclDataStart = D.GetexclDataStart(range); % find correct data, where to cut spectrum

%% Fit -> fit results are in property: A.FitResult
D.Fit;

% fit again with averaged response function
%D.LoadSingleRunObj;
%%
%     RF = zeros(D.nRuns,D.ModelObj.nTe,D.ModelObj.nqU);
%     for i=1:A.nRuns
%         RF(i,:,:) = D.SingleRunObj{i}.RF(1:D.ModelObj.nTe,:);
%     end
% %RFmean = D.SingleRunObj.
%% Display result
D.PlotFit('saveplot','OFF','FitResultsFlag','OFF');

