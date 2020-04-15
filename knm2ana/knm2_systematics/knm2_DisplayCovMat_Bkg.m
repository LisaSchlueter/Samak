% Plot Covariance Matrix for background (qU) sloep systematic effect
% systematics setting
RecomputeFlag = 'OFF';
nTrials       = 1000; 
SysBudget     = 34; % 31= knm2 preliminary input

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
    'TwinBias_Q',18573.7,...  % twin endpoint
    'SysBudget',SysBudget,...
    'RingMerge','Full',...
    'NonPoissonScaleFactor',1};  



CmArg = {'BkgCM','ON',...
    'SysEffects',struct('BkgShape','ON'),...
    'RecomputeFlag',RecomputeFlag,...
    'nTrials',nTrials,...
    };
%% init model
M = MultiRunAnalysis(RunArg{:});
M.exclDataStart = M.GetexclDataStart(40);
Sigma = std(E0);
FSDArg = {'SanityPlot','OFF','Sigma',Sigma};
M.ModelObj.LoadFSD(FSDArg{:});
M.ModelObj.ComputeTBDDS; M.ModelObj.ComputeTBDIS;

M.Fit;
FitResultsStat = M.FitResult;
%% calculate covariance matrix
M.chi2 = 'chi2CMShape';
M.ComputeCM(CmArg{:});

SysErrdef = GetSysErr(SysBudget);
M.InitModelObj_Norm_BKG;
M.FitCM_Obj.ComputeCM_Background('Display','ON','MaxSlopeCpsPereV',5.2*1e-06,'Mode','Gauss');

M.Fit;
FitResultsCM = M.FitResult;
%% display results 
mNuSysNeg =  sqrt(FitResultsCM.errNeg(1)^2-FitResultsStat.errNeg(1)^2);
mNuSysPos =  sqrt(FitResultsCM.errPos(1)^2-FitResultsStat.errPos(1)^2);
mNuSys =  0.5*(mNuSysNeg+mNuSysPos);
fprintf('mnuSq sensitivity stat only        = %.3f eV^2 \n',0.5*(-FitResultsStat.errNeg(1)+FitResultsStat.errPos(1)));
fprintf('mnuSq sensitivity stat + syst only = %.3f eV^2 \n',0.5*(-FitResultsCM.errNeg(1)+FitResultsCM.errPos(1)));
fprintf('mnuSq sensitivity syst only        = %.3f eV^2 \n',mNuSys);

%% plot cov mat
%M.FitCM_Obj.PlotCM('qUWindowIndexMax',-20,'saveplot','ON');




