% Calculte and/or plot systematics breakdown for KNM2
% MultiRing!
% MC figure skating
range = 40;

E0 = knm2FS_GetE0Twins('SanityPlot','OFF','Mode','FS2');
chi2 = 'chi2Stat';

%% set up model
RunArg = {'RunList','KNM2_Prompt',...
    'chi2','chi2Stat',...
    'DataType','Twin',...
    'fixPar','mNu E0 Bkg Norm',...
    'RadiativeFlag','ON',...
    'minuitOpt','min ; minos',...
    'FSDFlag','BlindingKNM2',...
    'ELossFlag','KatrinT2A20',...
    'SysBudget',37,...
    'AnaFlag','Ring',...
    'RingMerge','Full',...
    'chi2','chi2Stat',...
    'pullFlag',4,...
    'TwinBias_Q',E0,...
    'ROIFlag','Default',...
    'MosCorrFlag','OFF',...
    'NonPoissonScaleFactor',1};

MR = MultiRunAnalysis(RunArg{:});
MR.exclDataStart = MR.GetexclDataStart(range);
if ~strcmp(chi2,'chi2Stat')
    MR.NonPoissonScaleFactor = 1.112;
    MR.SetNPfactor; % convert to right dimension (if multiring)
    MR.chi2 = chi2;
   % MR.ComputeCM;
end

%%  broadening of FSD
if numel(E0)>0
    FSDArg = {'SanityPlot','OFF','Sigma',std(E0)};
    MR.ModelObj.LoadFSD(FSDArg{:});
end
%%
SysAll    = {'TASR','LongPlasma','Bkg'}; %Bkg has to be last
SysLeg    = {'Tritium activity fluctuations';'Long. source potential';'Background slope'};
S = RunSensitivity('RunAnaObj',MR,'SysEffectsAll',SysAll,'SysEffectLeg',SysLeg);
S.RecomputeFlag='OFF';
S.LimitFlag = 'Central';
S.ConfLevel=0; % 0 == 1 sigma
%%
S.PlotSysBreakdownBars2('Ranges',MR.exclDataStart,'SavePlot','pdf','HoldOn','OFF','SysInfoBox','OFF');

