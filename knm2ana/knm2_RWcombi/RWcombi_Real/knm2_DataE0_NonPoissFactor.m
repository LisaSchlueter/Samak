range = 40; % fit range in eV below E0
RunAnaArg = {'fixPar','E0 Bkg Norm',...         % free Parameter !!
    'DataType','Real',...              % Real, Twin or Fake
    'FSDFlag','BlindingKNM2',...       % final state distribution (theoretical calculation)
    'ELossFlag','KatrinT2',...         % energy loss function     ( different parametrizations available)
    'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
    'chi2','chi2Stat',...              % statistics only
    'RingMerge','Full'};

meanE0Err  = zeros(3,1);
dNormSigma = zeros(3,1);
%% RW1
M1 = MultiRunAnalysis('RunList','KNM2_RW1',RunAnaArg{:});
M1.exclDataStart = M1.GetexclDataStart(range);
M1.FitRunList;
E01    = M1.SingleRun_FitResults.chi2Stat.E0;
E0Err1 = M1.SingleRun_FitResults.chi2Stat.E0Err;
meanE01 = wmean(E01,E0Err1);
meanE0Err(1) =  mean(E0Err1);
dNorm1 = fitdist(E01-meanE01,'Normal');
dNormSigma(1) = dNorm1.sigma;
NP1 = dNorm1.sigma/meanE0Err(1);

%% RW2
M2 = MultiRunAnalysis('RunList','KNM2_RW2',RunAnaArg{:});
M2.exclDataStart = M2.GetexclDataStart(range);
M2.FitRunList;
E02    = M2.SingleRun_FitResults.chi2Stat.E0;
E0Err2 = M2.SingleRun_FitResults.chi2Stat.E0Err;
meanE02 = wmean(E02,E0Err2);
meanE0Err(2) =  mean(E0Err2);
dNorm2 = fitdist(E02-meanE02,'Normal');
dNormSigma(2) = dNorm2.sigma;
NP2 = dNorm2.sigma/meanE0Err(2);

%% RW3
M3 = MultiRunAnalysis('RunList','KNM2_RW3',RunAnaArg{:});
M3.exclDataStart = M3.GetexclDataStart(range);
M3.FitRunList;
E03    = M3.SingleRun_FitResults.chi2Stat.E0;
E0Err3 = M3.SingleRun_FitResults.chi2Stat.E0Err;
meanE03 = wmean(E03,E0Err3);
meanE0Err(3) =  mean(E0Err3);
dNorm3 = fitdist(E03-meanE03,'Normal');
dNormSigma(3) = dNorm3.sigma;
NP3 = dNorm3.sigma/meanE0Err(3);

%%

