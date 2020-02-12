% twins for the three periods seperatetly 
RunList = 'KNM2_Prompt';
fixPar = 'mNu E0 Bkg Norm'; % free parameter
DataType = 'Twin';
FSDFlag = 'BlindingKNM2';
ELossFlag = 'KatrinT2';
AnaFlag = 'StackPixel'; % uniform FPD
exclDataStart = 11;

chi2 = 'chi2Stat';
RunAnaArg = {'RunList',RunList,...
    'fixPar',fixPar,...
    'DataType',DataType,...
    'FSDFlag',FSDFlag,...
    'ELossFlag',ELossFlag,...
    'AnaFlag',AnaFlag,'chi2',chi2,...
    'exclDataStart',exclDataStart};

% read data and set up model
T = MultiRunAnalysis(RunAnaArg{:},'TwinBias_Q','Fit');
T.Fit;

TE = MultiRunAnalysis(RunAnaArg{:},'TwinBias_Q',18573.7);



