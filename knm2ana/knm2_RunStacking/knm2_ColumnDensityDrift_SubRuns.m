% calculate Twim simulation with actual subrun column density drift
% fit with average response function
% uniform FPD
range     = 40;
freePar   = 'mNu E0 Bkg Norm';
chi2      = 'chi2Stat';
DataType  = 'Twin';
AnaFlag   = 'StackPixel';
RingMerge = 'None';
DopplerEffectFlag = 'FSD';
BKG_PtSlope = 3*1e-06;
TwinBias_BKG_PtSlope = 3*1e-06;
FSDFlag   = 'KNM2_0p1eV';
SysBudget = 40;

savedir = [getenv('SamakPath'),'knm2ana/knm2_RunStack/results/'];
savename = sprintf('%sknm2_ColumnDensityDrift_SubRuns_Bpng-%.1fmucpsPers_%s_%.0feV_%s_%s_%s_%s.mat',...
    savedir,BKG_PtSlope*1e6,DataType,range,strrep(freePar,' ',''),chi2,AnaFlag,FSDFlag);

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

if exist(savename,'file')
    load(savename,'FitResult','RunAnaArg','A');
else
    SigmaSq =  0.0124+0.0025;
 
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
        'BKG_PtSlope',BKG_PtSlope,...
        'TwinBias_BKG_PtSlope',TwinBias_BKG_PtSlope,...
        'DopplerEffectFlag',DopplerEffectFlag};
    A = MultiRunAnalysis(RunAnaArg{:});
    %%
    A.exclDataStart = A.GetexclDataStart(range);
    nqUFit = numel(A.RunData.qU(A.exclDataStart:end));
    WGTS_CD_MolPerCm2_i = A.RunData.WGTS_CD_MolPerCm2; % average column density
    WGTS_CD_MolPerCm2_SubRun = A.RunData.WGTS_CD_MolPerCm2_SubRun(A.exclDataStart:end);
    
    RF_average      = A.ModelObj.RF;
    RF_subrunCD = zeros(size(RF_average,1),size(RF_average,2),nqUFit);
    
    % calculate RF for every subrun column density
    for i=1:nqUFit
        A.ModelObj.WGTS_CD_MolPerCm2 = WGTS_CD_MolPerCm2_SubRun(i);
        A.ModelObj.InitializeRF;
        RF_subrunCD(:,:,i) = A.ModelObj.RF;
    end
    
    %% construct new response function
    RF_subrun = zeros(size(RF_average));
      for i=1:nqUFit
         RF_subrun(:,i+A.exclDataStart-1) = squeeze(RF_subrunCD(:,i+A.exclDataStart-1,i));
      end
    %%
    A.ModelObj.RF = RF_subrun;
    A.ModelObj.ComputeTBDDS('B_bias',0.2216-A.ModelObj.BKG_RateSec_i);
    A.ModelObj.ComputeTBDIS;
    TBDIS_subrun = A.ModelObj.TBDIS; 
    A.RunData.TBDIS = TBDIS_subrun;
    A.Fit;
    FitResult_subrun = A.FitResult;
    
     A.ModelObj.RF = RF_average;
     A.Fit;
     FitResult_average = A.FitResult;
   
    MakeDir(savedir);
    save(savename,...
    'WGTS_CD_MolPerCm2_i','WGTS_CD_MolPerCm2_SubRun',...
    'RF_subrun','RF_average',...
    'FitResult_subrun','FitResult_average','RunAnaArg')
end
%%



