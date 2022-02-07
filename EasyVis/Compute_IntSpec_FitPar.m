% knm2: calculate spectral model (uniform) used for final knm2 fit

nIteration = 1000;
savedir = [getenv('SamakPath'),'EasyVis/results/'];
savename = sprintf('%sComputeIntSpec_FitPar_%.0f.mat',savedir,nIteration);


if exist(savename,'file') 
    load(savename,'FitResult','RunAnaArg','A');
else
    SigmaSq =  0.0124+0.0025;
    range     = 40;
    freePar   = 'mNu E0 Bkg Norm';
    chi2      = 'chi2CMShape';
    NonPoissonScaleFactor = 1.112;
    DataType  = 'Real';
    AnaFlag   = 'StackPixel';
    RingMerge = 'None';
    DopplerEffectFlag = 'FSD';
    BKG_PtSlope = 3*1e-06;
    TwinBias_BKG_PtSlope = 3*1e-06;
    FSDFlag   = 'KNM2_0p1eV';
    PullFlag = 99;%[7,24]; %24 = 3.0 mucps/s
    SysBudget = 40;
    AnaStr = AnaFlag;
    
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
    FitResult = A.FitResult;
    %%
    qU = A.ModelObj.qU(A.exclDataStart:end);
    Time = A.ModelObj.qUfrac(A.exclDataStart:end).*A.ModelObj.TimeSec;
    TBDIS_bf = A.ModelObj.TBDIS(A.exclDataStart:end);
    TBDIS_data = A.RunData.TBDIS(A.exclDataStart:end);
    
    MakeDir(savedir);
    save(savename,'FitResult','RunAnaArg','A','SigmaSq',...
        'Time','qU','TBDIS_bf','TBDIS_data')
    
    %% illustrate fit parameter: signal normalization
    N  = FitResult.par(4)+linspace(-0.7,0.7,nIteration);
    TBDIS_N = zeros(numel(TBDIS_bf),numel(N));
    for i=1:numel(N)
        A.ModelObj.ComputeTBDDS(...
            'mSq_bias',FitResult.par(1),...
            'E0_bias',FitResult.par(2),...
            'B_bias',FitResult.par(3),...
            'N_bias',N(i))
        A.ModelObj.ComputeTBDIS;
        TBDIS_N(:,i) = A.ModelObj.TBDIS(A.exclDataStart:end);
    end
    
%     % sanity plot
%     GetFigure;
%     plot(qU,TBDIS_N./Time);
%     set(gca,'YScale','log');
%     hold on;
%     plot(qU,TBDIS_bf./Time,'-k','LineWidth',2);
%     
    %% illustrate fit parameter: background
    B  = FitResult.par(3)+linspace(-0.050,0.050,nIteration);
    TBDIS_B = zeros(numel(TBDIS_bf),numel(N));
    for i=1:numel(B)
        progressbar(i/nIteration)
        A.ModelObj.SetFitBias(0);
        A.ModelObj.ComputeTBDDS(...
            'mSq_bias',FitResult.par(1),...
            'E0_bias',FitResult.par(2),...
            'B_bias',B(i),...
            'N_bias',FitResult.par(4))
        A.ModelObj.ComputeTBDIS;
        TBDIS_B(:,i) = A.ModelObj.TBDIS(A.exclDataStart:end);
    end
    
    % sanity plot
%     GetFigure;
%     plot(qU,TBDIS_B./Time);
%     set(gca,'YScale','log');
%     hold on;
%     plot(qU,TBDIS_bf./Time,'-k','LineWidth',2);
    % illustrate fit parameter: endpoint
    E0  = FitResult.par(2)+linspace(-5,5,nIteration);
    TBDIS_E0 = zeros(numel(TBDIS_bf),numel(N));
    for i=1:numel(B)
        A.ModelObj.SetFitBias(0);
        A.ModelObj.ComputeTBDDS(...
            'mSq_bias',FitResult.par(1),...
            'E0_bias',E0(i),...
            'B_bias',FitResult.par(3),...
            'N_bias',FitResult.par(4))
        A.ModelObj.ComputeTBDIS;
        TBDIS_E0(:,i) = A.ModelObj.TBDIS(A.exclDataStart:end);
    end
    
    % sanity plot
    GetFigure;
%     plot(qU,TBDIS_E0./Time);
%     set(gca,'YScale','log');
%     hold on;
%     plot(qU,TBDIS_bf./Time,'-k','LineWidth',2);
    % illustrate fit parameter: nu-mass
    mNuSq  = FitResult.par(2)+linspace(-50,50,nIteration);
    TBDIS_mNuSq = zeros(numel(TBDIS_bf),numel(N));
    for i=1:numel(mNuSq)
        A.ModelObj.SetFitBias(0);
        A.ModelObj.ComputeTBDDS(...
            'mSq_bias',mNuSq(i),...
            'E0_bias',FitResult.par(2),...
            'B_bias',FitResult.par(3),...
            'N_bias',FitResult.par(4))
        A.ModelObj.ComputeTBDIS;
        TBDIS_mNuSq(:,i) = A.ModelObj.TBDIS(A.exclDataStart:end);
    end
    
    % sanity plot
%     GetFigure;
%     plot(qU,TBDIS_mNuSq./Time);
%     set(gca,'YScale','log');
%     hold on;
%     plot(qU,TBDIS_bf./Time,'-k','LineWidth',2);
     
    save(savename,'FitResult','RunAnaArg','A','SigmaSq',...
        'Time','qU','TBDIS_bf','TBDIS_data',...
        'TBDIS_N','TBDIS_B','TBDIS_E0','TBDIS_mNuSq',...
        'N','B','E0','mNuSq');
end
%%