% knm2: calculate spectral model (uniform) used for final knm2 fit
% draw some random fit parameter values and calculate chi2

savedir = [getenv('SamakPath'),'EasyVis/results/'];
savename = sprintf('%sComputeIntSpec_Chi2.mat',savedir);


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
    TBDIS_data_err = sqrt(TBDIS_data);
    
    MakeDir(savedir);
    save(savename,'FitResult','RunAnaArg','A','SigmaSq',...
        'Time','qU','TBDIS_bf','TBDIS_data')
    
    %% illustrate fit parameter: signal normalization
    N  = FitResult.par(4)+linspace(-0.7,0.7,100);
    B  = FitResult.par(3)+linspace(-0.050,0.050,100);
    E0  = FitResult.par(2)+linspace(-5,5,100);
    mNuSq  = FitResult.par(2)+linspace(-50,50,100);
    
    %% combine randomly
    N = N(randi(100,100,1));
    B = B(randi(100,100,1));
    E0 = E0(randi(100,100),1);
    mNuSq = mNuSq(randi(100,100),1);
    %
          
    TBDIS = zeros(numel(TBDIS_bf),numel(N)); 
    chi2min = zeros(numel(N),1);

    for i=1:numel(N)
        A.ModelObj.ComputeTBDDS(...
            'mSq_bias',mNuSq(i),...
            'E0_bias',E0(i),...
            'B_bias',B(i),...
            'N_bias',N(i));
        A.ModelObj.ComputeTBDIS;
        TBDIS(:,i) = A.ModelObj.TBDIS(A.exclDataStart:end);
        chi2min(i) = sum((TBDIS(:,i)-TBDIS_data).^2./TBDIS_data);
    end
    
    chi2_bf = sum((TBDIS_bf-TBDIS_data).^2./TBDIS_data);
    
    % sort: ascending chi2
    [chi2min, Idx] = sort(chi2min);
   TBDIS = TBDIS(:,Idx);
    N = N(Idx);
    B = B(Idx);
    E0 = E0(Idx);
    mNuSq = mNuSq(Idx);
    
    %%
    % sanity plot
    GetFigure;
    plot(qU,TBDIS./Time);
    set(gca,'YScale','log');
    hold on;
    plot(qU,TBDIS_bf./Time,'-k','LineWidth',2);
    
    
    save(savename,'FitResult','RunAnaArg','A','SigmaSq',...
        'Time','qU','TBDIS_bf','TBDIS_data',...
        'TBDIS','N','B','E0','mNuSq','chi2min','chi2_bf');
end
%%