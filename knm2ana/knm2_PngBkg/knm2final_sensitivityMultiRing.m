% Calculte and/or plot systematics breakdown for KNM2 unblinding III
range = 40;
DataType = 'Real';
SigmaSq =  0.0124+0.0025;
PullFlag = 99;
freePar = 'mNu E0 Bkg Norm qU';
if contains(freePar,'qU')
    PullFlag = 99;
end
DopplerEffectFlag = 'FSD';
BKG_PtSlope = 3*1e-06;
TwinBias_BKG_PtSlope = 3*1e-06;
SysBudget = 41;

RunAnaArg = {'RunList','KNM2_Prompt',... % all KNM2 golden runs
    'fixPar',freePar,...           % free Parameter !!
    'DataType',DataType,...
    'FSDFlag','KNM2',...       % final state distribution (theoretical calculation)
    'ELossFlag','KatrinT2A20',...         % energy loss function     ( different parametrizations available)
    'AnaFlag','Ring',...
    'RingMerge','Full',...
    'chi2','chi2Stat',...              % statistics only
    'NonPoissonScaleFactor',1.112,...
    'TwinBias_Q',18573.70,...
    'SysBudget',SysBudget,...
    'chi2','chi2Stat',...
    'TwinBias_FSDSigma',sqrt(SigmaSq),...
    'FSD_Sigma',sqrt(SigmaSq),...
    'PullFlag',PullFlag,...
    'DopplerEffectFlag',DopplerEffectFlag,...
    'BKG_PtSlope',BKG_PtSlope,...
    'TwinBias_BKG_PtSlope',TwinBias_BKG_PtSlope};

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
