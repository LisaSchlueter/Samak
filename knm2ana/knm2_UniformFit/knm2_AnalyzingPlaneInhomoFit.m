% Get std of Bana and qU for kNM2 golden runs

savedir = [getenv('SamakPath'),'knm2ana/knm2_UniformFit/results/'];
savename = sprintf('%sknm2_AnalyzingPlaneInhomoFit.mat',savedir);

if exist(savename,'file')
    load(savename)
else
    MakeDir(savedir);
    
    range     = 40;
    freePar   = 'mNu E0 Bkg Norm';
    AnaFlag   = 'StackPixel';
    DopplerEffectFlag = 'FSD';
    BKG_PtSlope = 3*1e-06;
    TwinBias_BKG_PtSlope = 3*1e-06;
    FSDFlag   = 'KNM2_0p1eV';
    
    SigmaSq =  0.0124+0.0025;
    RunAnaArg = {'RunList','KNM2_Prompt',...
        'DataType','Real',...
        'fixPar',freePar,...
        'RadiativeFlag','ON',...
        'minuitOpt','min ; minos',...
        'FSDFlag',FSDFlag,...
        'ELossFlag','KatrinT2A20',...
        'AnaFlag',AnaFlag,...
        'chi2','chi2Stat',...
        'TwinBias_Q',18573.7,...
        'NonPoissonScaleFactor',1,...
        'FSD_Sigma',sqrt(SigmaSq),...
        'TwinBias_FSDSigma',sqrt(SigmaSq),...
        'PullFlag',99,...;%99 = no pull
        'BKG_PtSlope',BKG_PtSlope,...
        'TwinBias_BKG_PtSlope',TwinBias_BKG_PtSlope,...
        'DopplerEffectFlag',DopplerEffectFlag};
    
    A           = MultiRunAnalysis(RunAnaArg{:});
    %%
    PixList     = A.PixList;
    nPix        = numel(PixList);
    
    RF_pixel    = cell(nPix,1);
    
    TBDIS_pixel = zeros(A.ModelObj.nqU,nPix);
    Ba_pixel  = zeros(1,nPix);
    qU_pixel  = zeros(A.ModelObj.nqU,nPix);
    Te_pixel  = cell(nPix,1);
    
    for i=1:nPix
        A.PixList = PixList(i);
        A.StackRuns; % only 1 pixel
        A.SimulateStackRuns;
        A.ModelObj.ComputeTBDDS('B_bias',0.2216-A.ModelObj.BKG_RateSec_i);
        A.ModelObj.ComputeTBDIS;
        
        TBDIS_pixel(:,i) = A.ModelObj.TBDIS;
        RF_pixel{i} = A.ModelObj.RF;
        qU_pixel(:,i) = A.ModelObj.qU;
        Ba_pixel(i) = A.ModelObj.MACE_Ba_T;
        Te_pixel{i} = A.ModelObj.Te;
    end
    
    save(savename,'TBDIS_pixel','RF_pixel','qU_pixel','Ba_pixel','Te_pixel');
    
    A           = MultiRunAnalysis(RunAnaArg{:}); % average model
    
    TBDIS = sum(TBDIS_pixel);
    A.RunData.TBDIS = TBDIS;
    A.Fit;
    Fit_average = A.FitResult;
    
    save(savename,'Fit_average','-append');
end  
