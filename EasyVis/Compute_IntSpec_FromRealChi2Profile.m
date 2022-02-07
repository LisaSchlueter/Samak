% knm2: calculate spectral model (uniform) used for final knm2 fit
% load chi^2 profile over nu-mass
% get nuisance parameter values for each point and calculate int. spec

savedir = [getenv('SamakPath'),'EasyVis/results/'];
savename = sprintf('%sComputeIntSpec_FromRealChi2Profile.mat',savedir);


if exist(savename,'file')
  %  load(savename,'FitResult','RunAnaArg','A');
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
    chi2_bf = FitResult.chi2min;
   
    qU = A.ModelObj.qU(A.exclDataStart:end);
    Time = A.ModelObj.qUfrac(A.exclDataStart:end).*A.ModelObj.TimeSec;
    TBDIS_bf = A.ModelObj.TBDIS(A.exclDataStart:end);
    TBDIS_data = A.RunData.TBDIS(A.exclDataStart:end);
    TBDIS_data_err = sqrt(TBDIS_data);
    %% load chi2 profile
    savenameChi2Profile = [getenv('SamakPath'),'tritium-data/fit/Knm2/Chi2Profile/Uniform/Chi2Profile_Real_UniformScan_mNu_Knm2_UniformFPD_chi2CMShape_SysBudget40_NP1.112_FitParE0BkgNorm_nFit50_KNM2_0p1eV_min-2.5_max2.5.mat'];
    d = importdata(savenameChi2Profile);

    %% load nuisance parameter
    % sort for ascending neutrino mass
    mNuSq  = reshape(d.ScanResults.ParScan,100,1);
    mNuSq(1) = []; % double
    [mNuSq , SortIdx]=sort(mNuSq,'ascend');
    
    E0  = reshape(squeeze(d.ScanResults.par(2,:,:)),100,1);
    E0(1) = [];% double
    E0 = E0(SortIdx);
    
    B  = reshape(squeeze(d.ScanResults.par(3,:,:)),100,1);
    B(1) = [];% double
    B = B(SortIdx);
    
    N  = reshape(squeeze(d.ScanResults.par(4,:,:)),100,1);
    N(1) = [];% double
    N = N(SortIdx);
    
    chi2min  = reshape(squeeze(d.ScanResults.chi2min),100,1);
    chi2min(1) = [];% double
    chi2min = chi2min(SortIdx);

   %%       
    TBDIS = zeros(numel(TBDIS_bf),numel(N)); 

    for i=1:numel(N)
        A.ModelObj.ComputeTBDDS(...
            'mSq_bias',mNuSq(i),...
            'E0_bias',E0(i),...
            'B_bias',B(i),...
            'N_bias',N(i));
        A.ModelObj.ComputeTBDIS;
        TBDIS(:,i) = A.ModelObj.TBDIS(A.exclDataStart:end);
    end
    
    
    save(savename,'FitResult','RunAnaArg','A','SigmaSq',...
        'Time','qU','TBDIS_bf','TBDIS_data',...
        'TBDIS','N','B','E0','mNuSq','chi2min','chi2_bf');
end
%%