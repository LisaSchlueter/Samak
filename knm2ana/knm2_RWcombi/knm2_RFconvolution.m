% test convolution of response function to account for Delta(qU) when stacking runs
% not pursued further
% use broadened FSDs instead

RunList = 'KNM2_RW1';
fixPar = 'mNu E0 Bkg Norm'; % free parameter
DataType = 'Twin';
FSDFlag = 'BlindingKNM2';
ELossFlag = 'KatrinT2';
AnaFlag = 'StackPixel'; % uniform FPD
exclDataStart = 14;
TwinBias_Q = 18573.7; %constant endpoint (absolute number doesn't matter here)
chi2 = 'chi2Stat';

CommonArg = {'RunList',RunList,...
    'fixPar',fixPar,...
    'DataType',DataType,...
    'FSDFlag',FSDFlag,...
    'ELossFlag',ELossFlag,...
    'AnaFlag',AnaFlag,'chi2',chi2,...
    'TwinBias_Q','Fit',...
    'exclDataStart',exclDataStart,...
    'TwinBias_Q',TwinBias_Q};
M = MultiRunAnalysis(CommonArg{:});
RF = M.ModelObj.RF; % response function size = nTe x nqU
%%
qU = M.RunData.qU; % central values
DeltaqU = std(M.SingleRunData.qU,0,2); %standard deviation of qU per subrun

gtmp = gaussian(M.ModelObj.Te,qU(10),DeltaqU(10));
g = (gtmp./simpsons(M.ModelObj.Te,gtmp));

RFconv = conv(RF(:,10),g);



