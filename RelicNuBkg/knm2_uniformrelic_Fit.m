% unblinded fit with penning trap background slope
range     = 40;
freePar   = 'mNu E0 Bkg Norm eta';   %qU for multiring fit
chi2      = 'chi2CMShape';
DataType  = 'Real';
AnaFlag   = 'StackPixel';   
RingMerge = 'None';%'None';     
DopplerEffectFlag = 'FSD';      %!different from KNM1
BKG_PtSlope = 3*1e-06;
TwinBias_BKG_PtSlope = 3*1e-06;
FSDFlag   = 'KNM2';     %'KNM2_0p5eV' for faster fit (rebinned FSDs)
PullFlag = 99;%[7,24]; %24 = 3.0 mucps/s; 99 for no pull
SysBudget = 40;

savedir = [getenv('SamakPath'),'knm2ana/knm2_PngBkg/results/'];
savename = sprintf('%sknm2relic_Fit_Bpng-%.1fmucpsPers_%s_%.0feV_%s_%s_%s.mat',...
    savedir,BKG_PtSlope*1e6,DataType,range,strrep(freePar,' ',''),chi2,FSDFlag);

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

if ~exist(savename,'file')
    load(savename,'FitResult','RunAnaArg','A');
else
    SigmaSq =  0.0124+0.0025;   %broadening durch longitudinale plasmainhomogenität + plasma drift
    
     if strcmp(RingMerge,'None') && strcmp(chi2,'chi2CMShape') && strcmp(AnaFlag,'Ring')
         chi2tmp = 'chi2Stat';
     else
         chi2tmp  = chi2;
     end
     
    RunAnaArg = {'RunList','KNM2_Prompt',...
        'DataType',DataType,...
        'fixPar',freePar,...
        'RadiativeFlag','ON',...
        'fitter','matlab',...
        'minuitOpt','min ; imp',...
        'FSDFlag',FSDFlag,...
        'ELossFlag','KatrinT2A20',...   %!different
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
        'DopplerEffectFlag',DopplerEffectFlag,...
        'AngularTFFlag','ON'};                  %KNM1 default: OFF
    A = MultiRunAnalysis(RunAnaArg{:});
    %%
    A.exclDataStart = A.GetexclDataStart(range);
    %A.exclDataStop  = A.ModelObj.nqU-1;
    
    if strcmp(DataType,'Twin')
        A.ModelObj.RFBinStep = 0.01;
        A.ModelObj.InitializeRF;
    end
    
    if  contains(freePar,'BkgPTSlope') && contains(freePar,'BkgSlope')  && strcmp(chi2,'chi2CMShape')
        A.ComputeCM('BkgPtCM','OFF','BkgCM','OFF');
    elseif  contains(freePar,'BkgSlope') && strcmp(chi2,'chi2CMShape')
         A.ComputeCM('BkgPtCM','ON','BkgCM','OFF');
    elseif contains(freePar,'BkgPTSlope') && strcmp(chi2,'chi2CMShape')
         A.ComputeCM('BkgPtCM','OFF','BkgCM','ON');
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
                        'BkgPtCM','ON',...
                        'BkgCM','ON');
    else
        
    end
    A.ModelObj.eta_i = -5.8e10;
    A.ModelObj.ComputeTBDDS;
    A.ModelObj.ComputeTBDIS;
    A.Fit;
    FitResult = A.FitResult;
    MakeDir(savedir);
    %save(savename,'FitResult','RunAnaArg','A','SigmaSq')
end
%%

%A.PlotFit;
fprintf('m_nu^2 = %.3f + %.3f %.3f eV^2       , ',FitResult.par(1),FitResult.errPos(1),FitResult.errNeg(1))
fprintf('mean err = %.3f eV^2 \n',(FitResult.errPos(1)-FitResult.errNeg(1))/2)
fprintf('E_0 = %.3f + %.3f eV  \n',FitResult.par(2)+A.ModelObj.Q_i,FitResult.err(2))
fprintf('chi2 = %.3f (%.0f dof), p = %.3f  \n',FitResult.chi2min,FitResult.dof,1-chi2cdf(FitResult.chi2min,FitResult.dof));
fprintf('eta = %.3fe10 + %.3fe10 - %.3fe10       , ',FitResult.par(18),abs(FitResult.errPos(5)),abs(FitResult.errNeg(5)))
fprintf('mean err = %.3fe10 \n',(abs(FitResult.errPos(5))+abs(FitResult.errNeg(5)))/2)

%%