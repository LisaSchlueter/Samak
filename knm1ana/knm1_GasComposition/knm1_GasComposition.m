

% compute covariance matrix for varyying FSD, based on trueness of LARA
% concentrations

% get correlation matrix
chi2 = 'chi2CMShape';
fitter = 'minuit';
savedir  = [getenv('SamakPath'),'knm1ana/knm1_GasComposition/results/'];
savefileResult = sprintf('%sknm1_GasComposition_%s_%s.mat',savedir,chi2,fitter);
if exist(savefileResult,'file')  && 1==2
%     dR = importdata(savefileResult);
load(savefileResult)
else 
    % Init Model Object and covariance matrix object
    R = MultiRunAnalysis('RunList','KNM1',...
        'chi2','chi2Stat',...
        'DataType','Twin',...
        'fixPar','mNu E0 Norm Bkg',... free parameter
        'RadiativeFlag','ON',...
        'NonPoissonScaleFactor',1,...
        'minuitOpt','min ; minos',...
        'fitter',fitter,...
        'FSDFlag','Sibille0p5eV',...
        'ELossFlag','KatrinT2',...
        'SysBudget',22,...
        'AngularTFFlag','OFF');
    R.exclDataStart = R.GetexclDataStart(40);  % 40eV range = 27 subruns
    savefileCovMat = sprintf('%sknm1_LARADataCovMat_5000samples.mat',savedir);
    dCov = importdata(savefileCovMat);
    
    R.RunData.TBDIS = dCov.TBDIS_i;%+R.ModelObj.BKG_RateSec.*R.ModelObj.qUfrac.*R.ModelObj.TimeSec;
    R.InitModelObj_Norm_BKG;
    R.Fit;
    
    FitResult_stat = R.FitResult;
    
    FitCMFrac_stat = R.FitCMFrac;
    FitCM_stat     = R.FitCM;
    
    % get covariance matrix
   % reset entries above endpoint
   qUIdx = find(R.ModelObj.qU-18573.7>0);
   dCov.covmatfrac(qUIdx,qUIdx) = 0;
   dCov.covmat(qUIdx,qUIdx) = 0;
   dCov.covmatfracshape(qUIdx,qUIdx) = 0;
   dCov.covmatshape(qUIdx,qUIdx) = 0;
  
   
    FitCMFrac = FitCMFrac_stat+ dCov.covmatfrac;
    FitCM = FitCM_stat+ dCov.covmat;
    FitCMFracShape = FitCMFrac_stat+ dCov.covmatfracshape;
    FitCMShape = FitCM_stat+ dCov.covmatshape;  
    
    R.FitCMFracShape = FitCMFracShape;
    R.FitCMShape = FitCMShape;
    R.FitCM = FitCM;
    R.FitCMFrac = FitCMFrac;
    
    R.chi2 = chi2;
 %   R.InitFitPar;
    R.Fit;
    FitResult_cm = R.FitResult;
    
    if strcmp(fitter,'minuit')
        SysErrNeg = sqrt(FitResult_cm.errNeg(1).^2-FitResult_stat.errNeg(1)^2);
        SysErrPos = sqrt(FitResult_cm.errPos(1).^2-FitResult_stat.errPos(1)^2);
        SysErrMean = 0.5.*(SysErrNeg+SysErrPos);
    elseif strcmp(fitter,'matlab')
        SysErrNeg = 0;
        SysErrPos = 0;
        SysErrMean = 0;
    end
    SysErr = sqrt(FitResult_cm.err(1).^2-FitResult_stat.err(1)^2);
    
    save(savefileResult,'FitResult_stat','FitResult_cm','FitResult_stat','FitCMFrac',...
        'SysErr','SysErrNeg','SysErrPos','SysErrMean');
end


