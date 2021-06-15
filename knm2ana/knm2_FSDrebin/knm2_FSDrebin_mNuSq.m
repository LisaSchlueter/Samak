% unblinded fit with penning track background slope
BinningFactor = 1;
FSDFlag   = 'KNM2_0p1eV_cut50eV';

%%
range     = 40;
freePar   = 'mNu E0 Bkg Norm';
chi2      = 'chi2Stat';
DataType  = 'Twin';
AnaFlag   = 'StackPixel';
RingMerge = 'Full';
DopplerEffectFlag = 'FSD';
BKG_PtSlope = 3*1e-06;
TwinBias_BKG_PtSlope = 3*1e-06;
PullFlag = 99;
SysBudget = 40;

savedir = [getenv('SamakPath'),'knm2ana/knm2_FSDrebin/results/'];
savename = sprintf('%sknm2_FSDrebin_Fit_%s_%.0feV_%s_%s_%s_%sx%.0f.mat',...
    savedir,DataType,range,strrep(freePar,' ',''),chi2,AnaFlag,FSDFlag,BinningFactor);

if strcmp(chi2,'chi2Stat')
    NonPoissonScaleFactor = 1;
elseif  strcmp(chi2,'chi2CMShape')
    NonPoissonScaleFactor = 1.112;
elseif strcmp(chi2,'chi2Stat+')
    NonPoissonScaleFactor = 1.112;
    chi2 = 'chi2Stat';
end
if ~strcmp(chi2,'chi2Stat')
    savename = strrep(savename,'.mat',sprintf('_SysBudget%.0f.mat',SysBudget));
end

if exist(savename,'file')  && 1==2
    load(savename,'FitResult','tFitSec','FSD_E','FSD_P');
else
    SigmaSq =  0.0124+0.0025;
    
     if strcmp(RingMerge,'None') && strcmp(chi2,'chi2CMShape')
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
     
     A.ModelObj.LoadFSD('BinningFactor',BinningFactor);
     FSD_E = A.ModelObj.TTexE;
     FSD_P = A.ModelObj.TTexP;
     
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
    
    tStart =tic;
    A.Fit;
    tFitSec = toc(tStart);
    FitResult = A.FitResult;
    MakeDir(savedir);
    save(savename,'FitResult','RunAnaArg','A','SigmaSq','tFitSec','BinningFactor','FSDFlag','FSD_E','FSD_P')
end
fprintf('-----------------------------------------------\n');
fprintf('FSDFlag=%s (%.0f bins), cputime=%.1fs \n',FSDFlag,numel(FSD_E),tFitSec)
fprintf('m_nu^2 = %.3g + %.3f %.3f eV^2       , ',FitResult.par(1),FitResult.errPos(1),FitResult.errNeg(1))
fprintf('mean err = %.3f eV^2 \n',(FitResult.errPos(1)-FitResult.errNeg(1))/2)
fprintf('E_0 = %.3f + %.3f eV  \n',FitResult.par(2)+A.ModelObj.Q_i,FitResult.err(2))
fprintf('chi2 = %.3f (%.0f dof), p = %.3f  \n',FitResult.chi2min,FitResult.dof,1-chi2cdf(FitResult.chi2min,FitResult.dof));
fprintf('-----------------------------------------------\n');

