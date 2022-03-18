% unblinded fit with penning trap background slope
range     = 40;
freePar   = 'mNu E0 Bkg Norm qU';
chi2      = 'chi2Stat';%CMShape';
DataType  = 'Real';
AnaFlag   = 'Ring';
RingMerge = 'Full';
DopplerEffectFlag = 'FSD';
BKG_PtSlope = 3*1e-06;
TwinBias_BKG_PtSlope = 3*1e-06;
FSDFlag   = 'KNM2_0p1eV';
PullFlag = 99;%99;%6;%[7,24]; %24 = 3.0 mucps/s

if strcmp(AnaFlag,'Ring')
    if strcmp(RingMerge,'Full')
        AnaStr = AnaFlag;
        SysBudget = 41;
    else
        AnaStr = sprintf('Ring%s',RingMerge);
        SysBudget = 40;
    end
else
    SysBudget = 40;
    AnaStr = AnaFlag;
end

if strcmp(chi2,'chi2Stat')
    NonPoissonScaleFactor = 1;
elseif  strcmp(chi2,'chi2CMShape')
    NonPoissonScaleFactor = 1.112;
elseif strcmp(chi2,'chi2Stat+')
    NonPoissonScaleFactor = 1.112;
    chi2 = 'chi2Stat';
end


SigmaSq =  0.0124+0.0025;

if strcmp(RingMerge,'None') && strcmp(chi2,'chi2CMShape') && strcmp(AnaFlag,'Ring')
    chi2tmp = 'chi2Stat';
else
    chi2tmp  = chi2;
end

RunAnaArg = {'RunList','KNM2_Prompt',...
    'DataType',DataType,...
    'fixPar',freePar,...
    'RadiativeFlag','ON',...
    'minuitOpt','min ; minos',...
    'FSDFlag',FSDFlag,...
    'ELossFlag','KatrinT2A20',...
    'SysBudget',SysBudget,...
    'AnaFlag',AnaFlag,...
    'chi2',chi2tmp,...
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

if strcmp(DataType,'Twin')
    A.ModelObj.RFBinStep = 0.01;
    A.ModelObj.InitializeRF;
end

if  contains(freePar,'BkgPTSlope') && contains(freePar,'BkgSlope')  && strcmp(chi2,'chi2CMShape')
    A.ComputeCM('BkgPTCM','OFF','BkgCM','OFF');
elseif  contains(freePar,'BkgSlope') && strcmp(chi2,'chi2CMShape')
    A.ComputeCM('BkgPTCM','ON','BkgCM','OFF');
elseif contains(freePar,'BkgPTSlope') && strcmp(chi2,'chi2CMShape')
    A.ComputeCM('BkgPTCM','OFF','BkgCM','ON');
end

if strcmp(RingMerge,'None') && strcmp(chi2,'chi2CMShape') && strcmp(AnaFlag,'Ring')
    A.chi2 = chi2;
    A.ComputeCM('SysEffects',  struct(...
        'RF_EL','OFF',...   % Response Function(RF) EnergyLoss
        'RF_BF','OFF',...   % RF B-Fields
        'RF_RX','OFF',...   % Column Density, inel cross ection
        'FSD','ON',...
        'TASR','ON',...
        'TCoff_RAD','OFF',...
        'TCoff_OTHER','ON',...
        'DOPoff','OFF',...
        'Stack','OFF',...
        'FPDeff','ON',...
        'LongPlasma','ON'),...
        'BkgPTCM','ON',...
        'BkgCM','ON');
end

ProfileResult = A.ComputeChi2Profile('Parameter','mNu','nFit',30,...
    'ParMin',-0.5,'ParMax',1);;





