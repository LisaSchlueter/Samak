% investigate possible impact of new syst. effect: Tritium on rear wall
% tritium spectrum modifications
% 1) different response function (starting pos.)
% 2) Endpoint shifted
% 3) Signal normalization
% 4) No additional background
% 5) maybe different FSD
ScaleRW = 0.58*1e-03;%0.012/1.181;%;
E0ShifteV = 0;%+1.64; % eV
B_RW_Source = 1.58;%1.23;%2.52;%
freePar = 'mNu E0 Norm Bkg';
range = 40;
savedir = [getenv('SamakPath'),'ksn2ana/ksn2_Systematics/results/'];
savename = [savedir,sprintf('ksn2_RWsyst_Fit%s_%.0feV_RWE0shift%.3geV_ScaleRWRate%.3g.mat',strrep(freePar,' ',''),range,E0ShifteV,ScaleRW)];

if B_RW_Source~=2.52 
   savename =  strrep(savename,'.mat',sprintf('_%.3gWGTS_B_T.mat',B_RW_Source));
end

if exist(savename,'file') 
    load(savename)
    fprintf('Load spectra from file %s \n',savename)
else

    RunAnaArg = {'RunList','KNM2_Prompt',...
        'DataType','Twin',...
        'fixPar',freePar,...%free par
        'SysBudget',40,...
        'fitter','minuit',...
        'minuitOpt','min;migrad',...
        'RadiativeFlag','ON',...
        'FSDFlag','KNM2_0p1eV',...
        'ELossFlag','KatrinT2A20',...
        'AnaFlag','StackPixel',...
        'chi2','chi2Stat',...
        'NonPoissonScaleFactor',1,...
        'FSD_Sigma',sqrt(0.0124+0.0025),...
        'TwinBias_FSDSigma',sqrt(0.0124+0.0025),...
        'TwinBias_Q',18573.7,...
        'PullFlag',99,...;%99 = no pull
        'BKG_PtSlope',3*1e-06,...
        'TwinBias_BKG_PtSlope',3*1e-06,...
        'DopplerEffectFlag','FSD'};
     A = MultiRunAnalysis(RunAnaArg{:});
     A.exclDataStart = A.GetexclDataStart(range);
     % init model and save init settings
     FitResults                 = A.InitModelObj_Norm_BKG;
     A.ModelObj.BKG_RateSec_i   = A.ModelObj.BKG_RateSec_i+FitResults.par(3);
     A.ModelObj.NormFactorTBDDS = A.ModelObj.NormFactorTBDDS.*(1+FitResults.par(4));
     NormFactorTBDDS_i          = A.ModelObj.NormFactorTBDDS;
     Bkg_i                      = A.ModelObj.BKG_RateSec_i;
     Q_i                        = A.ModelObj.Q_i;
     WGTS_CD_MolPerCm2_i        = A.ModelObj.WGTS_CD_MolPerCm2;
    
     % compute TBDIS for regular tritium spectrum (wgts)
     A.ModelObj.ComputeTBDDS; A.ModelObj.ComputeTBDIS; 
     TBDIS_wgts             = A.ModelObj.TBDIS;
     %% init rear wall spectrum
     A.ModelObj.KTFFlag = 'RW_WGTSMACE';
     A.ModelObj.WGTS_B_T = B_RW_Source;
     A.ModelObj.ComputeNormFactorTBDDS;
     A.ModelObj.InitializeRF;
     
     NormFactorTBDDS_rw       = A.ModelObj.NormFactorTBDDS; % different norm factor, becauseof B-source
     A.ModelObj.Q_i           = Q_i + E0ShifteV;
     A.ModelObj.BKG_RateSec_i = 0; % no (additional) background
     A.ModelObj.BKG_PtSlope_i = 0;
     A.ModelObj.NormFactorTBDDS = NormFactorTBDDS_rw*ScaleRW; % scale to lower activity
     
     A.ModelObj.ComputeTBDDS; A.ModelObj.ComputeTBDIS;
     TBDIS_rw = A.ModelObj.TBDIS;
    
     %% reset model to regular settings
     A.ModelObj.BKG_PtSlope_i = 3*1e-06;
     A.ModelObj.BKG_RateSec_i = Bkg_i;
     A.ModelObj.NormFactorTBDDS = NormFactorTBDDS_i;
     A.ModelObj.Q_i = Q_i;
     A.ModelObj.KTFFlag = 'WGTSMACE';
     A.ModelObj.WGTS_B_T = 2.52;
     A.ModelObj.WGTS_CD_MolPerCm2 = WGTS_CD_MolPerCm2_i;
     A.ModelObj.ComputeNormFactorTBDDS;
     A.ModelObj.InitializeRF;
 
     % fit for wgts
     A.RunData.TBDIS = TBDIS_wgts;
     A.Fit;
     FitResults_ref = A.FitResult;
     
     % fit to wgts + rw 
     TBDIS_sum =   TBDIS_wgts+TBDIS_rw;
     A.RunData.TBDIS = TBDIS_sum;
     A.Fit;
     FitResults = A.FitResult;
     
     % save
     MakeDir(savedir);
     qU = A.RunData.qU;
     Time = A.RunData.qUfrac.*A.RunData.TimeSec;
     
     save(savename,'qU','TBDIS_sum','TBDIS_rw','TBDIS_wgts','E0ShifteV','Q_i','Time','ScaleRW','FitResults','FitResults_ref');
end
 %% display
fprintf('Fit results: \n')
fprintf('Ref fit: m^2 = %.3g eV^2 , chi2min = %.2g \n',FitResults_ref.par(1),FitResults_ref.chi2min);
fprintf('rw  fit: m^2 = %.3g eV^2 , chi2min = %.2g \n',FitResults.par(1),FitResults.chi2min);