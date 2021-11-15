

% compute covariance matrix for varyying FSD, based on trueness

% get correlation matrix
savedir  = [getenv('SamakPath'),'knm1ana/knm1_GasComposition/results/'];
savefileResult = sprintf('%sknm1_GasComposition.mat',savedir);
if exist(savefileResult,'file') & 1==2
    dR = importdata(savefileResult);
else
    
    % Init Model Object and covariance matrix object
    R = MultiRunAnalysis('RunList','KNM1',...
        'chi2','chi2Stat',...
        'DataType','Twin',...
        'fixPar','mNu E0 Norm Bkg',... free parameter
        'RadiativeFlag','ON',...
        'NonPoissonScaleFactor',1,...
        'minuitOpt','min ; minos',...
        'FSDFlag','SibilleFull',...
        'ELossFlag','KatrinT2',...
        'SysBudget',22,...
        'AngularTFFlag','OFF');
    R.exclDataStart = R.GetexclDataStart(40);  % 40eV range = 27 subruns
    R.Fit;
    FitResult_stat = R.FitResult;
    
    FitCMFrac_stat = R.FitCMFrac;
    FitCM_stat = R.FitCM;
    
    % get covariance matrix
    nSamples = 5e3;
    savefileCovMat = sprintf('%sknm1_LARADataCovMat_%.0fsamples.mat',savedir,nSamples);
    dCov = importdata(savefileCovMat);
    
    FitCMFrac = FitCMFrac_stat+ dCov.covmatfrac;
    FitCM = FitCM_stat+ dCov.covmat;
    R.chi2 = 'chi2CMShape';
    R.FitCMFracShape = FitCMFrac;
    R.FitCMShape = FitCM;
    R.Fit;
    FitResult_cm = R.FitResult;
    
    SysErrNeg = sqrt(FitResult_cm.errNeg(1).^2-FitResult_stat.errNeg(1)^2);
    SysErrPos = sqrt(FitResult_cm.errPos(1).^2-FitResult_stat.errPos(1)^2);
    SysErrMean = 0.5.*(SysErrNeg+SysErrPos);
    SysErr = sqrt(FitResult_cm.err(1).^2-FitResult_stat.err(1)^2);
    
    save(savefileResult,'FitResult_stat','FitResult_cm','FitResult_stat','FitCMFrac',...
        'SysErr','SysErrNeg','SysErrPos','SysErrMean');
end


