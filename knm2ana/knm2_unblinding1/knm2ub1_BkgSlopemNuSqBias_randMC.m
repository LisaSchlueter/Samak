%% randomized MC study
% stat + NP, 40 eV 
% 1. fit (stat + NP) with fixed bkg slope
% 2. fit (stat + NP + bkgslope cm) with fixed bkg slope
% compare results in scatter plot
range   = 40;
nSamples = 1e3;

savedir = [getenv('SamakPath'),'knm2ana/knm2_unblinding1/results/BestFit/'];
savename = sprintf('%sknm2ub1_BkgSlopeFreeVsFixmNuSqBiasRandMC_%.0feV_%.0fsamples.mat',savedir,range,nSamples);

if exist(savename,'file')
    load(savename);
else
    SigmaSq =  0.0124+0.0025;
    
    RunAnaArg = {'RunList','KNM2_Prompt',...
        'chi2','chi2Stat',...
        'DataType','Twin',...
        'fixPar','mNu E0 Bkg Norm',...
        'RadiativeFlag','ON',...
        'minuitOpt','min ; minos',...
        'FSDFlag','BlindingKNM2',...
        'ELossFlag','KatrinT2A20',...
        'chi2','chi2Stat',...
        'SysBudget',38,...
        'AnaFlag','StackPixel',...
        'TwinBias_Q',18573.7,...
        'NonPoissonScaleFactor',1.112,...
        'FSD_Sigma',sqrt(SigmaSq),...
        'TwinBias_FSDSigma',sqrt(SigmaSq)};
    
    %% Object:  stat + NP , bkg slope fixed to 0
    A = MultiRunAnalysis(RunAnaArg{:});
    A.exclDataStart = A.GetexclDataStart(range);
    A.ModelObj.RFBinStep = 0.01;
    A.ModelObj.InitializeRF;
    %% Object:  stat + NP + bkg slope cov mat
    B = MultiRunAnalysis(RunAnaArg{:});
    B.exclDataStart = B.GetexclDataStart(range);
    B.ModelObj.RFBinStep = 0.01;
    B.ModelObj.InitializeRF;
    B.chi2 = 'chi2CMShape';
    B.ComputeCM('BkgCM','ON','SysEffects',struct('FSD','OFF'));
    %% randomize integral spectrum
    TBDIS_i = A.RunData.TBDIS;
    [StatCM, StatCMFrac]         =  A.ComputeCM_StatPNP;
    A.FitCM     = StatCM;       A.FitCMShape = StatCM;
    A.FitCMFrac = StatCMFrac;   A.FitCMFracShape = StatCMFrac;
    TBDIS = mvnrnd(TBDIS_i',A.FitCM,nSamples)';
    
    %% make copy of object for parallel computing
    Apar = copy(repmat(A,nSamples,1));
    Bpar = copy(repmat(B,nSamples,1));
    
    %A.RunData.TBDIS = TBDIS_i+sqrt(TBDIS_i).*randn(A.ModelObj.nqU,1); % randomize
    mNuSq      = zeros(nSamples,1);
    mNuSqErr   = zeros(nSamples,1);
    FitResults = cell(nSamples,1);
    
    mNuSqCM      = zeros(nSamples,1);
    mNuSqErrCM   = zeros(nSamples,1);
    FitResultsCM = cell(nSamples,1);
  %%  
    parfor i=1:nSamples
        Apar(i).SimulateStackRuns;
        Apar(i).RunData.TBDIS = TBDIS(:,i);
        %% stat. only - fixed background slope
        Apar(i).Fit;
        mNuSq(i)      = Apar(i).FitResult.par(1);
        mNuSqErr(i)   = (Apar(i).FitResult.errPos(1)-Apar(i).FitResult.errNeg(1))/2;
        FitResults{i} = Apar(i).FitResult;
        %% stat. only - background slope cm
        Bpar(i).SimulateStackRuns;
        Bpar(i).RunData.TBDIS = TBDIS(:,i);
        Bpar(i).Fit;
        mNuSqCM(i)      = Bpar(i).FitResult.par(1);
        mNuSqErrCM(i)   = (Bpar(i).FitResult.errPos(1)-Bpar(i).FitResult.errNeg(1))/2;
        FitResultsCM{i} = Bpar(i).FitResult;
    end
    
    save(savename,'mNuSq','mNuSqErr','FitResults','mNuSqCM','mNuSqErrCM','FitResultsCM','RunAnaArg');
end

%% plot
% p = scatterhist(mNuSq,mNuSqCM,'Direction','out',...);%,'HistogramDisplayStyle','bar',...
%     'Location','Northeast','Color',rgb('DodgerBlue'));
GetFigure
h1 = histogram(mNuSq-mNuSqCM,'Normalization','probability',...
    'FaceColor',rgb('SkyBlue'),'FaceAlpha',1,'EdgeColor',rgb('PowderBlue'));
PrettyFigureFormat('FontSize',22);
xlabel(sprintf('\\Delta{\\itm}_\\nu^2 (eV^{ 2})'));
ylabel('Frequency')
title(sprintf('stat. only - stat. and {\\itB}^{slope} syst. , %.0f samples',nSamples),'FontWeight','normal','FontSize',get(gca,'FontSize')-2);
leg = legend(sprintf(' \\mu = %.3f eV , \\sigma = %.3f eV ',mean(mNuSq-mNuSqCM),std(mNuSq-mNuSqCM)));
leg.EdgeColor = rgb('Silver'); 
leg.Location='northwest';

xlim([-0.1 0.1]);
plotdir = [getenv('SamakPath'),'knm2ana/knm2_unblinding1/plots/'];
plotname = sprintf('%sknm2ub1_BkgSlopeFreeVsFixmNuSqBiasRandMC_%.0feV_%.0fsamples.png',plotdir,range,nSamples);
print(plotname,'-dpng','-r350');

