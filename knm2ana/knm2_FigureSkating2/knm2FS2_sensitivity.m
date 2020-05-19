% Calculte and/or plot systematics breakdown for KNM2 figure skating II

% MC figure skating
range = 40;

E0 = knm2FS_GetE0Twins('SanityPlot','OFF','Mode','FS2');

RunAnaArg = {'RunList','KNM2_Prompt',... % all KNM2 golden runs
    'fixPar','mNu E0 Bkg Norm',...           % free Parameter !!
    'DataType','Twin',...              
    'FSDFlag','BlindingKNM2',...       % final state distribution (theoretical calculation)
    'ELossFlag','KatrinT2A20',...         % energy loss function     ( different parametrizations available)
    'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
    'chi2','chi2Stat',...              % statistics only
    'NonPoissonScaleFactor',1.112,...
    'TwinBias_Q',E0,...
    'ROIFlag','Default',...
    'SysBudget',36,...
    'chi2','chi2Stat'};

MC = MultiRunAnalysis(RunAnaArg{:});
MC.exclDataStart = MC.GetexclDataStart(range);
%%
if numel(E0)>0
    FSDArg = {'SanityPlot','OFF','Sigma',std(E0)};
    MC.ModelObj.LoadFSD(FSDArg{:});
end
%%
SysAll    = {'TASR','LongPlasma','Bkg'}; %Bkg has to be last
SysLeg    = {'Tritium activity fluctuations';'Long. source potential';'Background slope'};

S = RunSensitivity('RunAnaObj',MC,'SysEffectsAll',SysAll,'SysEffectLeg',SysLeg);
S.RecomputeFlag='OFF';
S.LimitFlag = 'Central';
S.ConfLevel=0; % 0 == 1 sigma
%%
S.AsymErr = 'ON';
S.PlotSysBreakdownBars2('Ranges',MC.exclDataStart,'SavePlot','pdf','HoldOn','OFF','SysInfoBox','OFF');

