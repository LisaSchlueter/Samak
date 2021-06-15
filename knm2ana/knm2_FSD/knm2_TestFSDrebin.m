savedir = [getenv('SamakPath'),'knm2ana/knm2_FSDrebin/results/'];

FSD_SigmaSq = 0.0124+0.0025;
BKG_PtSlope = 3*1e-06;

RunAnaArg = {'RunList','KNM2_Prompt',...
    'DataType','Twin',...
    'fixPar','mNu E0 Norm Bkg',...
    'minuitOpt','min ; minos',...
    'ELossFlag','KatrinT2A20',...
    'SysBudget',40,...
    'AnaFlag','StackPixel',...
    'chi2','chi2Stat',...
    'TwinBias_Q',18573.7,...
    'NonPoissonScaleFactor',1,...
    'FSD_Sigma',sqrt(FSD_SigmaSq),...
    'TwinBias_FSDSigma',sqrt(FSD_SigmaSq),...
    'RingMerge','Full',...
    'PullFlag',99,...
    'BKG_PtSlope',BKG_PtSlope,...
    'TwinBias_BKG_PtSlope',BKG_PtSlope,...
    'DopplerEffectFlag','FSD'};
%% not rebinned FSDs
file1 = sprintf('%sknm2_TestFSDrebin_KNM2.mat',savedir);
if exist(file1,'file')
    d1 = importdata(file1);
    FitResult1 = d1.FitResult;
    nBins1     = d1.nBins;
    FitTime1   = d1.FitTime;
else
    T1 = MultiRunAnalysis(RunAnaArg{:},'FSDFlag','KNM2','exclDataStart',11);
    tStart = tic; 
    T1.Fit;
    FitTime = toc(tStart);
    FitResult = T1.FitResult;
    FSDFlag = T1.FSDFlag;
    
    FSD_Energy = T1.ModelObj.TTexE;
    FSD_Prob   = T1.ModelObj.TTexP;
    nBins      = numel(FSD_Energy);
    MakeDir(savedir);
    save(file1,'FitTime','FitResult','RunAnaArg','FSDFlag','FSD_Energy','FSD_Prob','nBins');
end
%% 0.1eV rebinned FSDs
file2 = sprintf('%sknm2_TestFSDrebin_KNM2_0p1eV.mat',savedir);
if exist(file2,'file')
    d2 = importdata(file2);
    FitResult2 = d2.FitResult;
    nBins2     = d2.nBins;
    FitTime2   = d2.FitTime;
else
    T2 = MultiRunAnalysis(RunAnaArg{:},'FSDFlag','KNM2_0p1eV','exclDataStart',11);
    tStart2 = tic;
    T2.Fit;
    FitTime = toc(tStart2);
    FitResult = T2.FitResult;
    FSDFlag   = T2.FSDFlag;
    
    FSD_Energy = T2.ModelObj.TTexE;
    FSD_Prob   = T2.ModelObj.TTexP;
    nBins      = numel(FSD_Energy);
    
    save(file2,'FitTime','FitResult','RunAnaArg','FSDFlag','FSD_Energy','FSD_Prob','nBins');
end
%% result
fprintf('Original FSDs (%.0f bins): mNuSq = %.4f eV^2 , fit time: %.1f s\n',nBins1,FitResult1.par(1),FitTime1);
fprintf('Rebinned FSDs (%.0f bins): mNuSq = %.4f eV^2 , fit time: %.1f s\n',nBins2,FitResult2.par(1),FitTime2);

%%
% T2 = RunAnalysis(RunAnaArg{:},'FSDFlag','KNM2_0p5eV','exclDataStart',11);
% tic; T2.Fit; toc;
% T3 = RunAnalysis(RunAnaArg{:},'FSDFlag','KNM2_0p1eV','exclDataStart',11);
% tic;  T3.Fit;toc;
    
    