% calculate twins with average endpoint from 3 RW setting periods
% runs from 1 period have common endpoint
% apply RW shift on model
% fist step: known shift from MC

% twins for the three periods seperatetly 
RunList = 'KNM2_Prompt';
FSDFlag = 'BlindingKNM2';
ELossFlag = 'KatrinT2';
AnaFlag = 'StackPixel'; % uniform FPD
exclDataStart = 11;

chi2 = 'chi2Stat';
RunAnaArg = {'FSDFlag',FSDFlag,...
    'ELossFlag',ELossFlag,...
    'AnaFlag',AnaFlag,'chi2',chi2,...
    'exclDataStart',exclDataStart};
%% Get stacked endpoint per RW period
D1 = MultiRunAnalysis(RunAnaArg{:},'RunList','KNM2_RW1','DataType','Real','fixPar','E0 Bkg Norm'); 
D1.Fit;
D2 = MultiRunAnalysis(RunAnaArg{:},'RunList','KNM2_RW2','DataType','Real','fixPar','E0 Bkg Norm'); 
D2.Fit;
D3 = MultiRunAnalysis(RunAnaArg{:},'RunList','KNM2_RW3','DataType','Real','fixPar','E0 Bkg Norm'); 
D3.Fit;
%% prepare twin bias
TwinBias_Q = [(D1.FitResult.par(2)+D1.ModelObj.Q_i).*ones(1,D1.nRuns),...
    (D2.FitResult.par(2)+D2.ModelObj.Q_i).*ones(1,D2.nRuns),...
    (D3.FitResult.par(2)+D3.ModelObj.Q_i).*ones(1,D3.nRuns,1)];
E0    = [D1.FitResult.par(2),D2.FitResult.par(2),D3.FitResult.par(2)]+D1.ModelObj.Q_i;
E0Err =[D1.FitResult.err(2),D2.FitResult.err(2),D3.FitResult.err(2)];
meanE0 = wmean(E0,1./E0Err.^2);

MultiGauss_Weights = [D1.nRuns,D2.nRuns,D3.nRuns]./sum([D1.nRuns,D2.nRuns,D3.nRuns]);
MultiGauss_RelPos  = E0-meanE0; %shift
%% calculte twins for three periods
T3P = MultiRunAnalysis(RunAnaArg{:},'RunList','KNM2_Prompt','TwinBias_Q',TwinBias_Q,'DataType','Twin','fixPar','mNu E0 Bkg Norm',...
    'Twin_SameqUFlag','ON'); 
%%
%T3P.ModelObj.LoadFSD('Sigma',0.01,'MultiPos',MultiGauss_RelPos,'MultiWeights',MultiGauss_Weights);%mean(std(T3P.SingleRunData.qU'))
T3P.Fit;


