% Calculte and/or plot systematics breakdown for KNM2
% MultiRing!
% MC figure skating
range = 40;

E0 = knm2FS_GetE0Twins('SanityPlot','OFF','Mode','FS2');
chi2 = 'chi2CMShape';

%% set up model
RunArg = {'RunList','KNM2_Prompt',...
    'chi2','chi2Stat',...
    'DataType','Twin',...
    'fixPar','mNu E0 Bkg Norm',...
    'RadiativeFlag','ON',...
    'minuitOpt','min ; minos',...
    'FSDFlag','BlindingKNM2',...
    'ELossFlag','KatrinT2A20',...
    'SysBudget',36,...
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
    MR.ComputeCM;
end

%%  broadening of RF + broadening/shift of FSD
TimeSec = zeros(3,1);
TimeSec(1) = sum(MR.SingleRunData.TimeSec(1:171));
TimeSec(2) = sum(MR.SingleRunData.TimeSec(172:268));
TimeSec(3) = sum(MR.SingleRunData.TimeSec(269:361));
MultiWeights = TimeSec./sum(TimeSec);
MultiPos = [E0(1),E0(end-120),E0(end)]';
MultiPosRel = repmat(MultiPos-wmean(MultiPos,MultiWeights),1,MR.nRings);
Sigma = repmat(std(E0),3,MR.nRings);
FSDArg = {'MultiPos',MultiPosRel,'MultiWeight',MultiWeights,...
    'SanityPlot','OFF','Sigma',Sigma};
MR.ModelObj.LoadFSD(FSDArg{:});

%%
S = RunSensitivity('RunAnaObj',MR);
S.RecomputeFlag='OFF';
S.LimitFlag = 'Central';
S.ConfLevel=0; % 0 == 1 sigma
%%
S.PlotSysBreakdownBars2('Ranges',MR.exclDataStart,'SavePlot','pdf','HoldOn','OFF','SysInfoBox','OFF');

