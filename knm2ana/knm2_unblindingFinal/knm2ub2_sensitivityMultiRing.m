% Calculte and/or plot systematics breakdown for KNM2 unblinding 1
range = 40;
DataType = 'Real';
SigmaSq =  0.0124+0.0025;
PullFlag = 99;
freePar = 'mNu E0 Bkg Norm qU';
if contains(freePar,'qU')
    PullFlag = 6;
end
RunAnaArg = {'RunList','KNM2_Prompt',... % all KNM2 golden runs
    'fixPar',freePar,...           % free Parameter !!
    'DataType',DataType,...
    'FSDFlag','Sibille0p5eV',...       % final state distribution (theoretical calculation)
    'ELossFlag','KatrinT2A20',...         % energy loss function     ( different parametrizations available)
    'AnaFlag','Ring',...
    'RingMerge','Full',...
    'chi2','chi2Stat',...              % statistics only
    'NonPoissonScaleFactor',1.112,...
    'TwinBias_Q',18573.70,...
    'SysBudget',41,...
    'chi2','chi2Stat',...
    'TwinBias_FSDSigma',sqrt(SigmaSq),...
    'FSD_Sigma',sqrt(SigmaSq),...
    'PullFlag',PullFlag};

MC = MultiRunAnalysis(RunAnaArg{:});
MC.exclDataStart = MC.GetexclDataStart(range);

%%
%SysAll    = {'TASR','FSD','Bkg'}; %Bkg has to be last
%SysLeg    = {'Tritium activity fluctuations';'FSD';'Background slope'};
S = RunSensitivity('RunAnaObj',MC);%,'SysEffectsAll',SysAll,'SysEffectLeg',SysLeg);
S.RecomputeFlag='OFF';
S.LimitFlag = 'Central';
S.ConfLevel=0; % 0 == 1 sigma
%%
S.AsymErr = 'ON';
S.PlotSysBreakdownBars2('Ranges',MC.exclDataStart,'SavePlot','pdf','HoldOn','OFF','SysInfoBox','OFF');
