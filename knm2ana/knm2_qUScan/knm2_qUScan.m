% Uniform fit on KNM2 data/twins
% January 2020, Lisa
ELossFlag = 'KatrinT2';%'Aseev'; 

RunAnaArg = {'RunList','KNM2_Prompt',...  % define run number -> see GetRunList
    'fixPar','mNu E0 Bkg Norm',...         % free Parameter !!
    'DataType','Twin',...              % Real, Twin or Fake
    'FSDFlag','BlindingKNM2',...       % final state distribution (theoretical calculation)
    'ELossFlag',ELossFlag,...         % energy loss function     ( different parametrizations available)
    'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
    'chi2','chi2Stat',...              % statistics only
    'NonPoissonScaleFactor',1,...
    'MosCorrFlag','OFF',...
    'TwinBias_Q',18573.7,...
    'ROIFlag','14keV',...
    'DopplerEffectFlag','FSD'};

%% build object of MultiRunAnalysis class
D = MultiRunAnalysis(RunAnaArg{:});

%% modify some parameters in your analysis
range = 40;               % fit range in eV below endpoint        
D.exclDataStart = D.GetexclDataStart(range); % find correct data, where to cut spectrum
%% extra label for qU-scan with "wrong" energy loss function
if ~strcmp(ELossFlag,'KatrinT2')
    saveStr = ['_',ELossFlag];
else
    saveStr = '';
end
%%
D.qUScan('qURange',[90 12],'RecomputeFlag','OFF',...
        'saveplot','ON','ErrorBarScaling',1,'saveStr',saveStr);
