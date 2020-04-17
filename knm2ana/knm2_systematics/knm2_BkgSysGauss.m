% background slope (qU) systematics
% 2 different strategies: fit & cut off or gaussian randomization
% systematics setting
RecomputeFlag = 'OFF';
CovMatRecomputeFlag = 'ON';
MaxSlopeCpsPereV = 13.3*1e-06;%5.2*1e-06;

savedir = [getenv('SamakPath'),'knm2ana/knm2_systematics/results/'];
MakeDir(savedir);
savename = sprintf('%sknm2_BKGsys_GaussOnly_%.1fmcpskeV.mat',savedir,MaxSlopeCpsPereV*1e6);
if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename)
else
% model setting
RunList   = 'KNM2_Prompt';
AnaFlag   = 'StackPixel'; % FPD segmentation
E0 = knm2FS_GetE0Twins('SanityPlot','OFF');

RunArg = {'FSDFlag','BlindingKNM2',...
    'ELossFlag','KatrinT2',...
    'AnaFlag',AnaFlag,...
    'chi2','chi2Stat',...      % switch on later by hand
    'RunList',RunList,...      
    'fixPar','mNu E0 Norm Bkg',... % free parameter
    'DataType','Twin',...      % MC twins
    'TwinBias_Q',18573.70,...  % twin endpoint
    'SysBudget',34,...
    'RingMerge','Full',...
    'NonPoissonScaleFactor',1};        

CmArg = {'BkgCM','ON',...
    'SysEffects',struct('BkgShape','ON'),...
    'nTrials',50000,...
    'MaxSlopeCpsPereV',MaxSlopeCpsPereV};
%% init model
M = MultiRunAnalysis(RunArg{:});
M.exclDataStart = M.GetexclDataStart(40);
Sigma = std(E0);
FSDArg = {'SanityPlot','OFF','Sigma',Sigma};
M.ModelObj.LoadFSD(FSDArg{:});
M.ModelObj.ComputeTBDDS; M.ModelObj.ComputeTBDIS;

%% stat. only
M.Fit;
FitResultsStat = M.FitResult;
%% stat. + syst: calculate covariance matrix: slope fit & cut off
M.chi2 = 'chi2CMShape';
%% stat. + syst: calculate covariance matrix: slope randn
M.ComputeCM(CmArg{:},'BkgMode','Gauss');
CMFracGauss      = M.FitCM_Obj.CovMatFrac;
CMFracGaussShape = M.FitCM_Obj.CovMatFracShape;

M.Fit;
FitResultsCMGauss = M.FitResult;

%% prepare results
mNuSysNegGauss =  sqrt(FitResultsCMGauss.errNeg(1)^2-FitResultsStat.errNeg(1)^2);
mNuSysPosGauss =  sqrt(FitResultsCMGauss.errPos(1)^2-FitResultsStat.errPos(1)^2);
mNuSysGauss =  0.5*(mNuSysNegGauss+mNuSysPosGauss);

save(savename,'FitResultsStat','FitResultsCMGauss','mNuSysGauss',...
               'CMFracGauss','CMFracGaussShape');
end
%% display results
mNuStatAsym = 0.5*(-FitResultsStat.errNeg(1)+FitResultsStat.errPos(1));
fprintf('mnuSq sensitivity stat only        = %.3f eV^2 \n',mNuStatAsym);
fprintf('mnuSq sensitivity stat + syst only = %.3f eV^2  (gauss) \n',0.5*(-FitResultsCMGauss.errNeg(1)+FitResultsCMGauss.errPos(1)));
fprintf(2,'mnuSq sensitivity syst only (asym) = %.3g eV^2  (gauss) \n',mNuSysGauss);
%% plot cov mat
%M.FitCM_Obj.PlotCM('qUWindowIndexMax',-20,'saveplot','ON');




