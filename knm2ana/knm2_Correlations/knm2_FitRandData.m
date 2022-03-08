%
% KNM2 Final fit Results
% Uniform Fit
% Golden Run List
% Golden Pixel List
% randomized twins

%% settings
nFit = 1;
RecomputeFlag = 'OFF';
range     = 40;
freePar   = 'mNu E0 Bkg Norm';
chi2      = 'chi2CMShape';
DataType  = 'Twin';
AnaFlag   = 'StackPixel';%'Ring';
RingMerge = 'None';
DopplerEffectFlag = 'FSD';
BKG_PtSlope = 3*1e-06;
TwinBias_BKG_PtSlope = 3*1e-06;
FSDFlag   = 'KNM2_0p1eV';
SysBudget    = 40;
NonPoissonScaleFactor = 1.112;
SigmaSq =  0.0124+0.0025;

savedir = [getenv('SamakPath'),'knm2ana/knm2_Correlations/results/'];
savefile = sprintf('%sknm2FitRandData_%s_%s_NP%.4g_%s_%.0feV_%s_%.0ffit.mat',...
    savedir,DataType,chi2,NonPoissonScaleFactor,strrep(freePar,' ',''),range,FSDFlag,nFit);

if strcmp(chi2,'chi2CMShape')
    savefile = strrep(savefile,chi2,sprintf('%s_SysBudget%.0f',chi2,SysBudget));
end

if exist(savefile,'file') && strcmp(RecomputeFlag,'OFF')
    load(savefile);
    fprintf('load from file %s \n',savefile)
else
    MakeDir(savedir);
    
     TwinmNuSq = 0.28;
     
    RunAnaArg = {'RunList','KNM2_Prompt',...
        'DataType',DataType,...
        'fixPar',freePar,...
        'RadiativeFlag','ON',...
        'DopplerEffectFlag',DopplerEffectFlag,...
        'minuitOpt','min ; minos',...
        'FSDFlag',FSDFlag,...
        'ELossFlag','KatrinT2A20',...
        'SysBudget',SysBudget,...
        'AnaFlag',AnaFlag,...
        'chi2',chi2,...
        'NonPoissonScaleFactor',NonPoissonScaleFactor,...
        'FSD_Sigma',sqrt(SigmaSq),...
        'RingMerge',RingMerge,...
        'PullFlag',99,...;%99 = no pull
        'BKG_PtSlope',BKG_PtSlope,...
        'TwinBias_Q',18573.7,...
        'TwinBias_FSDSigma',sqrt(SigmaSq),...
        'TwinBias_mnuSq',TwinmNuSq,...
        'TwinBias_BKG_PtSlope',TwinBias_BKG_PtSlope};
    
    A = MultiRunAnalysis(RunAnaArg{:});
    A.exclDataStart = A.GetexclDataStart(range);
    A.InitModelObj_Norm_BKG('RecomputeFlag','ON');
    A.ComputeCM;
   
    [FitPar, FitErr, FitChi2min, dof,TBDIS]  = A.FitTwin('nSamples',nFit);
    
    MakeDir(savedir);
    save(savefile,'A','FitPar', 'FitErr', 'FitChi2min', 'dof','TBDIS');
end
