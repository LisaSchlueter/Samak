% Calculte and/or plot systematics breakdown for KNM2 unblinding 1
range = 40;
DataType = 'Real';
SigmaSq =  0.0124+0.0025;
SysBudget = 40;
BKG_PtSlope = 3*1e-06;
TwinBias_BKG_PtSlope = 3*1e-06; 
DopplerEffectFlag = 'FSD';

RunAnaArg = {'RunList','KNM2_Prompt',... % all KNM2 golden runs
    'fixPar','mNu E0 Bkg Norm',...           % free Parameter !!
    'DataType',DataType,...
    'FSDFlag','KNM2',...       % final state distribution (theoretical calculation)
    'ELossFlag','KatrinT2A20',...         % energy loss function     ( different parametrizations available)
    'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
    'chi2','chi2Stat',...              % statistics only
    'NonPoissonScaleFactor',1.112,...
    'TwinBias_Q',18573.70,...
    'SysBudget',SysBudget,...
    'chi2','chi2Stat',...
    'TwinBias_FSDSigma',sqrt(SigmaSq),...
    'FSD_Sigma',sqrt(SigmaSq),...
    'BKG_PtSlope',BKG_PtSlope,...
    'TwinBias_BKG_PtSlope',TwinBias_BKG_PtSlope,...
    'DopplerEffectFlag',DopplerEffectFlag};

MC = MultiRunAnalysis(RunAnaArg{:});
MC.exclDataStart = MC.GetexclDataStart(range);

if strcmp(DataType,'Twin')
    MC.ModelObj.RFBinStep = 0.01;
    MC.ModelObj.InitializeRF;
end

%%
%SysAll    = {'TASR','FSD','Bkg'}; %Bkg has to be last
%SysLeg    = {'Tritium activity fluctuations';'FSD';'Background slope'};

S = RunSensitivity('RunAnaObj',MC);%,'SysEffectsAll',SysAll,'SysEffectLeg',SysLeg);
S.RecomputeFlag='OFF';
S.LimitFlag = 'Central';
S.ConfLevel=0; % 0 == 1 sigma
%%
AsymErr = 'ON';
S.PlotSysBreakdownBars2('Ranges',MC.exclDataStart,'SavePlot','pdf','HoldOn','OFF','SysInfoBox','OFF');
