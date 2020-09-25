SigmaSq =  0.0124+0.0025;

RunAnaArg = {'RunList','KNM2_Prompt',...
    'chi2','chi2Stat',...
    'DataType','Twin',...
    'fixPar',freePar,...
    'RadiativeFlag','ON',...
    'minuitOpt','min ; minos',...
    'FSDFlag','BlindingKNM2',...
    'ELossFlag','KatrinT2A20',...
    'SysBudget',SysBudget,...
    'AnaFlag','StackPixel',...
    'RingMerge','Full',...
    'chi2',chi2,...
    'pullFlag',99,...
    'TwinBias_Q',18573.7,...
    'NonPoissonScaleFactor',1.112,...
    'TwinBias_FSDSigma',sqrt(SigmaSq)};
A = MultiRunAnalysis(RunAnaArg{:});
A.exclDataStart = A.GetexclDataStart(range);