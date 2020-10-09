% Calculte and/or plot energy-loss systematics breakdown for KNM2

range = 40;
mNuSqErr = zeros(4,1);
RunAnaArg = {'RunList','KNM2_Prompt',... % all KNM2 golden runs
    'fixPar','mNu E0 Bkg Norm',...           % free Parameter !!
    'DataType','Twin',...
    'FSDFlag','BlindingKNM2',...       % final state distribution (theoretical calculatin
    'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
    'chi2','chi2Stat',...              % statistics only
    'NonPoissonScaleFactor',1,...
    'ROIFlag','Default',...
    'SysBudget',36,...
    'chi2','chi2Stat'};

M = MultiRunAnalysis(RunAnaArg{:},'ELossFlag','KatrinT2A20','TwinBias_Q',18573.56);
M.exclDataStart = M.GetexclDataStart(range);
TBDIS_i = M.RunData.TBDIS;
%% stat only
M.Fit;
mNuSqErrStat = M.FitResult.err(1);

%% stat + syst. KatrinT2A20
M.chi2 = 'chi2CMShape';
M.ComputeCM('SysEffects',struct('RF_EL','ON'),'BkgCM','OFF')
M.Fit;
mNuSqErr(1) = M.FitResult.err(1);

%% stat + syst. KatrinT2
M1 = MultiRunAnalysis(RunAnaArg{:},'ELossFlag','KatrinT2','TwinBias_Q',18573.55);
M1.exclDataStart = M.GetexclDataStart(range);
M1.chi2 = 'chi2CMShape';
M1.RunData.TBDIS = TBDIS_i;
M1.ComputeCM('SysEffects',struct('RF_EL','ON'),'BkgCM','OFF')
M1.Fit;
mNuSqErr(2) = M1.FitResult.err(1);

%% stat + syst. Aseev
M2 = MultiRunAnalysis(RunAnaArg{:},'ELossFlag','Aseev','TwinBias_Q',18573.54);
M2.exclDataStart = M2.GetexclDataStart(range);
M2.chi2 = 'chi2CMShape';
M2.RunData.TBDIS = TBDIS_i;
M2.ComputeCM('SysEffects',struct('RF_EL','ON'),'BkgCM','OFF')
M2.Fit;
mNuSqErr(3) = M2.FitResult.err(1);

%% stat + syst. Abdurashitov
M3 = MultiRunAnalysis(RunAnaArg{:},'ELossFlag','Abdurashitov','TwinBias_Q',18573.53);
M3.exclDataStart = M3.GetexclDataStart(range);
M3.chi2 = 'chi2CMShape';
M3.RunData.TBDIS = TBDIS_i;
M3.ComputeCM('SysEffects',struct('RF_EL','ON'),'BkgCM','OFF')
M3.Fit;
mNuSqErr(4) = M3.FitResult.err(1);

%% syst only
mNuSqSys = mNuSqErr.^2-mNuSqErrStat^2;

