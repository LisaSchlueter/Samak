% calculate chi^2 for KNM1 and KNM2 using best fit parameter of COMMON fit
DataType = 'Real';
chi2 = 'chi2CMShape';
mNuSqCommon = 0.10;
KNM1SysBudget = 24;
KNM1Doppler   = 'OFF';

% KNM1
savefile = [getenv('SamakPath'),...
    sprintf('knm2ana/knm2_Combination/results/knm2_Combination_Chi2CommonPar_Knm1SysBudget%.0f_Knm1DE%s_mNuSqCommon%.3geV2.mat',...
    KNM1SysBudget,KNM1Doppler,mNuSqCommon)];

if exist(savefile,'file') 
    load(savefile)
else
    %% knm1
    range   = 40;         % 40eV range = 27 subruns
    Real = MultiRunAnalysis('RunList','KNM1',...
        'chi2',chi2,...
        'DataType',DataType,...
        'fixPar','mNu E0 Norm Bkg',... free parameter
        'RadiativeFlag','ON',...
        'NonPoissonScaleFactor',1.064,...
        'minuitOpt','min ; minos',...
        'FSDFlag','Sibille0p5eV',...
        'ELossFlag','KatrinT2',...
        'SysBudget',KNM1SysBudget,...
        'AngularTFFlag','OFF',...
        'DopplerEffectFlag',KNM1Doppler);
    Real.exclDataStart = Real.GetexclDataStart(range);
    Real.Fit;
    TBDIS_DataKnm1 = Real.RunData.TBDIS(Real.exclDataStart:end);
    FitResult_Knm1 = Real.FitResult;
    
    % check
    Real.ModelObj.ComputeTBDDS(...
        'mSq_bias',Real.FitResult.par(1),...
        'E0_bias',Real.FitResult.par(2),...
        'B_bias',Real.FitResult.par(3),...
        'N_bias',Real.FitResult.par(4));
    Real.ModelObj.ComputeTBDIS;
    TBDIS_ModelKnm1     = Real.ModelObj.TBDIS(Real.exclDataStart:end);
    CovMatShape_Knm1 = Real.FitCMShape(Real.exclDataStart:end,Real.exclDataStart:end);
    chi2min_Knm1    = chi2func(TBDIS_DataKnm1,TBDIS_ModelKnm1,CovMatShape_Knm1);
    
    % calculate for common best fit par
    Real.fixPar = ConvertFixPar('freePar','E0 Bkg Norm','nPar',Real.nPar,'nPixel',Real.ModelObj.nPixels);
    Real.ModelObj.mnuSq_i = mNuSqCommon;
    Real.Fit;
    FitResult_Knm1Common = Real.FitResult;
    TBDIS_ModelKnm1Common     = Real.ModelObj.TBDIS(Real.exclDataStart:end);
    chi2min_Knm1Common   = chi2func(TBDIS_DataKnm1,TBDIS_ModelKnm1Common,CovMatShape_Knm1);
    
    save(savefile,'FitResult_Knm1','CovMatShape_Knm1','chi2min_Knm1','TBDIS_DataKnm1','TBDIS_ModelKnm1',...
        'TBDIS_ModelKnm1Common','chi2min_Knm1Common','FitResult_Knm1Common','mNuSqCommon');
    
    %% knm2
    freePar               = 'mNu E0 Bkg Norm';
    AnaFlag               = 'StackPixel';
    RingMerge             = 'Full';%'None'; %only relevand when AnaFlag = Ring
    DopplerEffectFlag     = 'FSD';
    BKG_PtSlope           = 3*1e-06;
    TwinBias_BKG_PtSlope  = 3*1e-06;
    FSDFlag               = 'KNM2_0p1eV';
    PullFlag              = 99;% 99 = no pull
    SysBudget             = 40;
    NonPoissonScaleFactor = 1.112;
    SigmaSq               =  0.0124+0.0025;
    
    % Init Model Object and covariance matrix object
    RunAnaArg = {'RunList','KNM2_Prompt',...
        'DataType',DataType,...
        'fixPar',freePar,...
        'RadiativeFlag','ON',...
        'minuitOpt','min ; minos',...
        'FSDFlag',FSDFlag,...
        'ELossFlag','KatrinT2A20',...
        'SysBudget',SysBudget,...
        'AnaFlag',AnaFlag,...
        'chi2',chi2,...
        'TwinBias_Q',18573.7,...
        'NonPoissonScaleFactor',NonPoissonScaleFactor,...
        'FSD_Sigma',sqrt(SigmaSq),...
        'TwinBias_FSDSigma',sqrt(SigmaSq),...
        'RingMerge',RingMerge,...
        'PullFlag',PullFlag,...;%99 = no pull
        'BKG_PtSlope',BKG_PtSlope,...
        'TwinBias_BKG_PtSlope',TwinBias_BKG_PtSlope,...
        'DopplerEffectFlag',DopplerEffectFlag};
    A = MultiRunAnalysis(RunAnaArg{:});
    A.exclDataStart = A.GetexclDataStart(range);
    
    A.Fit;
    TBDIS_DataKnm2 = A.RunData.TBDIS(A.exclDataStart:end);
    FitResult_Knm2 = A.FitResult;
    
    % check
    A.ModelObj.ComputeTBDDS(...
        'mSq_bias',A.FitResult.par(1),...
        'E0_bias',A.FitResult.par(2),...
        'B_bias',A.FitResult.par(3),...
        'N_bias',A.FitResult.par(4));
    A.ModelObj.ComputeTBDIS;
    TBDIS_ModelKnm2     = A.ModelObj.TBDIS(A.exclDataStart:end);
    CovMatShape_Knm2 = A.FitCMShape(A.exclDataStart:end,A.exclDataStart:end);
    chi2min_Knm2    = chi2func(TBDIS_DataKnm2,TBDIS_ModelKnm2,CovMatShape_Knm2);
    
    % calculate for common best fit par
    A.fixPar = ConvertFixPar('freePar','E0 Bkg Norm','nPar',A.nPar,'nPixel',A.ModelObj.nPixels);
    A.ModelObj.mnuSq_i = mNuSqCommon;
    A.Fit;
    FitResult_Knm2Common  = A.FitResult;
    TBDIS_ModelKnm2Common = A.ModelObj.TBDIS(A.exclDataStart:end);
    chi2min_Knm2Common    = chi2func(TBDIS_DataKnm2,TBDIS_ModelKnm2Common,CovMatShape_Knm2);
    
    save(savefile,'FitResult_Knm2','CovMatShape_Knm2','chi2min_Knm2','TBDIS_DataKnm2','TBDIS_ModelKnm2',...
        'TBDIS_ModelKnm2Common','chi2min_Knm2Common','FitResult_Knm2Common','mNuSqCommon','-append');
    
end

function chi2 = chi2func(y,m,CovMatShapeFit)
chi2 = ((y - m)')* (CovMatShapeFit  \ (y - m));
end