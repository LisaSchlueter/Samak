% syst. uncertainty contribution from background time slope (penning trap)
% unblinded fit with penning track background slope
PullTerm = 'ON';
PullFlag = 24; % only used in fit with free Bpt slope: 24==3e-6 cps/s
DataType  = 'Twin';
BKG_PtSlopeErr = 3*1e-06;
BKG_PtSlope = 3*1e-06;
TwinBias_BKG_PtSlope = 3*1e-06;
RunList = 'KNM2_Prompt';
range     = 40;
freePar   = 'mNu E0 Bkg Norm';
AnaFlag   = 'StackPixel';
RingMerge = 'Full';
FSDFlag   = 'KNM2';
SigmaSq =  0.0124+0.0025;

%%
savedir      = [getenv('SamakPath'),'knm2ana/knm2_PngBkg/results/'];
savenameCM   = sprintf('%sknm2_BkgPT_SysErr_%s_%s_BPT%.1fmuCpsPerS_BPTerr%.1fmuCpsPerS_CM.mat',savedir,DataType,strrep(freePar,' ',''),1e6*BKG_PtSlope,1e6*BKG_PtSlopeErr);
savenameStat = sprintf('%sknm2_BkgPT_SysErr_%s_%s_BPT%.1fmuCpsPerS_Stat.mat',savedir,DataType,strrep(freePar,' ',''),1e6*BKG_PtSlope);

if strcmp(DataType','Twin')
    savenameCM = strrep(savenameCM,'.mat',sprintf('_TwinBPT%.1fmuCpsPerS.mat',1e6*TwinBias_BKG_PtSlope));
    savenameStat = strrep(savenameStat,'.mat',sprintf('_TwinBPT%.1fmuCpsPerS.mat',1e6*TwinBias_BKG_PtSlope));
end

RunAnaArg = {'RunList','KNM2_Prompt',...
    'DataType',DataType,...
    'fixPar',freePar,...
    'RadiativeFlag','ON',...
    'minuitOpt','min ; minos',...
    'FSDFlag',FSDFlag,...
    'ELossFlag','KatrinT2A20',...
    'SysBudget',40,...
    'AnaFlag',AnaFlag,...
    'TwinBias_Q',18573.7,...
    'NonPoissonScaleFactor',1,...
    'FSD_Sigma',sqrt(SigmaSq),...
    'TwinBias_FSDSigma',sqrt(SigmaSq),...
    'RingMerge',RingMerge,...
    'PullFlag',99,...%99 = no pull
    'BKG_PtSlope',BKG_PtSlope,...
    'TwinBias_BKG_PtSlope',TwinBias_BKG_PtSlope,...
    'DopplerEffectFlag','FSD',...
    'chi2','chi2Stat'};


if exist(savenameStat,'file') && 1==2 
    fprintf('load %s \n',savenameStat);
    load(savenameStat,'mNuSqErrStat','FitResultStat','A');
else
    A = MultiRunAnalysis(RunAnaArg{:});
    A.exclDataStart = A.GetexclDataStart(range);
    A.Fit;
    FitResultStat = A.FitResult;
    mNuSqErrStat = 0.5*(A.FitResult.errPos(1)- A.FitResult.errNeg(1));
    save(savenameStat,'mNuSqErrStat','FitResultStat','RunAnaArg','A');
end


if exist(savenameCM,'file') 
    fprintf('load %s \n',savenameCM);
    load(savenameCM,'mNuSqErrCM','FitResultCM');
else
    A.chi2 = 'chi2CMShape';
    % compute ONLY background PT covariance matrix
    A.ComputeCM('BkgPTCM','ON','RecomputeFlag','ON','BKG_PtSlopeErr',BKG_PtSlopeErr,...
        'SysEffect',struct('FSD','OFF'),'BkgCM','OFF','nTrials',1e4);
    A.Fit;
    mNuSqErrCM = 0.5*(A.FitResult.errPos(1)- A.FitResult.errNeg(1));
    FitResultCM = A.FitResult;
    save(savenameCM,'mNuSqErrCM','FitResultCM','BKG_PtSlopeErr','RunAnaArg');
    
    % reset
    A.chi2 = 'chi2Stat';
end


if strcmp(PullTerm,'ON') && any(PullFlag==[24,25])
    savenamePull = strrep(savenameCM,'CM.mat',sprintf('Pull%.0f.mat',PullFlag));
    if exist(savenamePull,'file')
        fprintf('load %s \n',savenamePull);
        load(savenamePull,'mNuSqErrPull','FitResultPull');
    else
        A.chi2 = 'chi2Stat';
        A.fixPar = ConvertFixPar('nPar',A.nPar,'freePar',[freePar,' BkgPTSlope']);
        A.pullFlag = PullFlag;
        A.Fit;
        
        mNuSqErrPull = 0.5*(A.FitResult.errPos(1)- A.FitResult.errNeg(1));
        FitResultPull = A.FitResult;
        
        save(savenamePull,'mNuSqErrPull','FitResultPull','BKG_PtSlopeErr','RunAnaArg');
        
        % reset pull & fit parameter
        A.pullFlag = 99;
        A.fixPar   = ConvertFixPar('nPar',A.nPar,'freePar',freePar);
    end
end

%% display
fprintf('slope uncertainty (%s ,%.1f mucps/s) 1 sigma sys err on m^2: \n',DataType,1e6.*BKG_PtSlopeErr);
fprintf('%.3f eV^2 (cov. mat) \n',sqrt(mNuSqErrCM^2-mNuSqErrStat^2));
fprintf('%.3f eV^2 (pull-term) \n',sqrt(mNuSqErrPull^2-mNuSqErrStat^2));
 