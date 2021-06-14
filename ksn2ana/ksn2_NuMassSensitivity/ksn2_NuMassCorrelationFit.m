% estimate correlation m4-m_nu
% method:
% 1. select MC truth for m_4,Ue4
% 2. vary m_4 slightly around MC truth & fit m_nu
DataType = 'Twin';
savedir = [getenv('SamakPath'),'ksn2ana/ksn2_NuMassSensitivity/results/CorrFits/'];
MakeDir(savedir);
nStepsAll = 2;
mNu4Sq_all = repmat(logspace(-1,log10(40^2),nStepsAll),nStepsAll,1);
mNu4Sq_all = reshape(mNu4Sq_all,nStepsAll^2,1);
sin2T4_all = repmat(logspace(-3,log10(0.5),nStepsAll)',1,nStepsAll);
sin2T4_all = reshape(sin2T4_all,nStepsAll^2,1);

for j=1:numel(mNu4Sq_all)
    mNu4Sq = mNu4Sq_all(j);
    sin2T4 = sin2T4_all(j);
    nSteps = 5;
    savefile = sprintf('%sksn2_NuMassCorrelationFit_%s_MCtruth_m4Sq%.3geV2_sin2T4%.3g_nSteps%.0f.mat',savedir,DataType,mNu4Sq,sin2T4,nSteps);
    if exist(savefile,'file')
        load(savefile);
    else
        %% settings that might change
        chi2 = 'chi2CMShape';
        nGridSteps = 50;
        range = 40;
        
        %% configure RunAnalysis object
        if strcmp(chi2,'chi2Stat')
            NonPoissonScaleFactor = 1;
        elseif  strcmp(chi2,'chi2CMShape')
            NonPoissonScaleFactor = 1.112;
        end
        RunAnaArg = {'RunList','KNM2_Prompt',...
            'DataType',DataType,...
            'fixPar','mNu E0 Norm Bkg',...%free par
            'SysBudget',40,...
            'fitter','minuit',...
            'minuitOpt','min;migrad',...
            'RadiativeFlag','ON',...
            'FSDFlag','KNM2_0p1eV',...
            'ELossFlag','KatrinT2A20',...
            'AnaFlag','StackPixel',...
            'chi2',chi2,...
            'NonPoissonScaleFactor',NonPoissonScaleFactor,...
            'FSD_Sigma',sqrt(0.0124+0.0025),...
            'TwinBias_FSDSigma',sqrt(0.0124+0.0025),...
            'TwinBias_Q',18573.7,...
            'PullFlag',99,...;%99 = no pull
            'BKG_PtSlope',3*1e-06,...
            'TwinBias_BKG_PtSlope',3*1e-06,...
            'DopplerEffectFlag','FSD'};
        A = MultiRunAnalysis(RunAnaArg{:});
        A.InitModelObj_Norm_BKG('RecomputeFlag','ON');
        %%
        
        A.ModelObj.SetFitBiasSterile(mNu4Sq,sin2T4);
        A.ModelObj.ComputeTBDDS(...
            'E0_bias',A.FitResult.par(2),...
            'B_bias',A.FitResult.par(3),...
            'N_bias',A.FitResult.par(4));
        A.ModelObj.ComputeTBDIS;
        TBDIS_i = A.ModelObj.TBDIS;
        A.RunData.TBDIS = TBDIS_i;
        
        
        mNu4Sq_test = mNu4Sq+linspace(-1,1,nSteps);
        FitResult = cell(nSteps,1);
        mNuSq = zeros(nSteps,1);
        for i=1:nSteps
            progressbar((i-1)./nSteps)
            A.ModelObj.SetFitBias(0);
            A.ModelObj.SetFitBiasSterile(mNu4Sq_test(i),sin2T4);
            A.Fit;
            FitResult{i} = A.FitResult;
            mNuSq(i) = A.FitResult.par(1);
        end
        
        CorrMat = corrcoef(mNu4Sq_test,mNuSq);
        CovMat  = cov(mNu4Sq_test,mNuSq);
        
        save(savefile,'mNuSq','mNu4Sq_test','mNu4Sq','sin2T4','FitResult','CorrMat','CovMat');
    end
end
