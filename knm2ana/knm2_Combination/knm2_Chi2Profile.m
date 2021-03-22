% KNM-2 chi2-profile

%% settings
range                 = 40;
freePar               = 'mNu E0 Bkg Norm qU';
chi2                  = 'chi2CMShape';
DataType              = 'Real';
AnaFlag               = 'Ring';
RingMerge             = 'Full';%'None'; %only relevand when AnaFlag = Ring
DopplerEffectFlag     = 'FSD';
BKG_PtSlope           = 3*1e-06;
TwinBias_BKG_PtSlope  = 3*1e-06;
FSDFlag               = 'KNM2';
PullFlag              = 99;% 99 = no pull
SysBudget             = 41;
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
%%
A.ResetFitResults;
A.FitResult.par(1) = -0.5;
ScanResults = A.GetAsymFitError('Mode','Uniform','ParScanMax',1.5,'nFitMax',20);
%%
mNuSq =[flipud(ScanResults.ParScan(:,2));ScanResults.ParScan(2:end,1)];
Chi2 = [flipud(ScanResults.chi2min(:,2));ScanResults.chi2min(2:end,1)];

plot(mNuSq,Chi2,':');