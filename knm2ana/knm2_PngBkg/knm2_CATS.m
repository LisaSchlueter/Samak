% unblinded fit with penning trap background slope
% CATS diagnostics

range     = 40;
freePar   = 'mNu E0 Bkg Norm';% qU';
chi2      = 'chi2Stat+';%CMShape';
DataType  = 'Real';
AnaFlag   = 'StackPixel';%Ring';
RingMerge = 'Full';
DopplerEffectFlag = 'FSD';
BKG_PtSlope = 3*1e-06;
TwinBias_BKG_PtSlope = 3*1e-06;
FSDFlag   = 'KNM2';

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

savedir = [getenv('SamakPath'),'knm2ana/knm2_PngBkg/results/'];
savename = sprintf('%sknm2ubfinal_Fit_Bpng-%.1fmucpsPers_%s_%.0feV_%s_%s_%s_%s_CATS.mat',...
    savedir,BKG_PtSlope*1e6,DataType,range,strrep(freePar,' ',''),chi2,AnaStr,FSDFlag);

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

if strcmp(DataType,'Twin')
    savename = strrep(savename,'.mat',sprintf('_TwinBpng-%.1fmucpsPers.mat',1e6*TwinBias_BKG_PtSlope));
end

if any(PullFlag~=99)
    if numel(PullFlag)==2
        savename = strrep(savename,'.mat',sprintf('_pull%.0f_%.0f.mat',PullFlag(1),PullFlag(2)));
    else
        savename = strrep(savename,'.mat',sprintf('_pull%.0f.mat',PullFlag));
    end
end

if exist(savename,'file') 
    load(savename,'FitResult','RunAnaArg','A','F');
else
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
    %%
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
    else
        
    end
    F = A.Fit('CATS','ON');
    FitResult = A.FitResult;
    MakeDir(savedir);
    save(savename,'FitResult','RunAnaArg','A','SigmaSq','F')
end

if ~exist('F','var')
    F =  A.Fit;
    save(savename,'F','-append');
end

%% fit result
fprintf('============================================\n');
fprintf('%s (%s) %s \n',A.AnaFlag,A.RingMerge,A.chi2)
fprintf('m_nu^2 = %.3f + %.2f %.2f eV^2, ',FitResult.par(1),FitResult.errPos(1),FitResult.errNeg(1))
fprintf('mean err = %.2f eV^2 \n',(FitResult.errPos(1)-FitResult.errNeg(1))/2)
fprintf('E_0 = %.3f + %.3f eV  \n',FitResult.par(2)+A.ModelObj.Q_i,FitResult.err(2))

fprintf('chi2 = %.3f (%.0f dof), p = %.3f  \n',FitResult.chi2min,FitResult.dof,1-chi2cdf(FitResult.chi2min,FitResult.dof));
fprintf('============================================\n ');

%%
%F.Samakcats_StandResidualLeverage2D('savePlot','OFF');
    F.Samakcats_Dfbetas_Abs('savePlot','OFF');








