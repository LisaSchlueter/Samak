% investigate possible impact of new syst. effect: Tritium on rear wall
% tritium spectrum modifications
% 1) different response function (starting pos.)
% 2) Endpoint shifted
% 3) Signal normalization
% 4) No additional background
% 5) maybe different FSD
ScaleRW = 1e-03;%5e-04;
E0ShifteV = 1.5; % eV
freePar = 'mNu E0 Norm Bkg';
savedir = [getenv('SamakPath'),'ksn2ana/ksn2_Systematics/results/'];
savename = [savedir,sprintf('ksn2_RWsyst_Fit%s_RWE0shift%.3geV_ScaleRWRate%.2g.mat',strrep(freePar,' ',''),E0ShifteV,ScaleRW)];

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
     FitResults= A.InitModelObj_Norm_BKG;
     A.ModelObj.BKG_RateSec_i = A.ModelObj.BKG_RateSec_i+FitResults.par(3);
     A.ModelObj.NormFactorTBDDS = A.ModelObj.NormFactorTBDDS.*(1+FitResults.par(4));
     A.ModelObj.ComputeTBDDS; A.ModelObj.ComputeTBDIS;
     
     TBDIS_wgts = A.ModelObj.TBDIS;
     Q_i = A.ModelObj.Q_i;
     
     
     %% init rear wall spectrum
     A.ModelObj.KTFFlag = 'RW_WGTSMACE';
     A.ModelObj.InitializeRF;
     A.ModelObj.Q_i = Q_i+E0ShifteV;
     Bkg_i = A.ModelObj.BKG_RateSec_i;
     A.ModelObj.BKG_RateSec_i = 0; % no (additional) background
     A.ModelObj.BKG_PtSlope_i = 0;
     A.ModelObj.ComputeTBDDS; A.ModelObj.ComputeTBDIS;
   
     TBDIS_rw = A.ModelObj.TBDIS;
     
     %% reset model to regular settings
     A.ModelObj.BKG_PtSlope_i = 3*1e-06;
     A.ModelObj.BKG_RateSec_i = Bkg_i;
     A.ModelObj.Q_i = Q_i;
     A.ModelObj.KTFFlag = 'WGTSMACE';
     A.ModelObj.InitializeRF;

     % fit for wgts
     A.RunData.TBDIS = TBDIS_wgts;
     A.Fit;
     FitResults_ref = A.FitResult;
     
     % fit to wgts + rw 
     TBDIS_sum =   TBDIS_wgts+TBDIS_rw.*ScaleRW;
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
fprintf('Ref fit: m^2 = %.2g eV^2 , chi2min = %.2g \n',FitResults_ref.par(1),FitResults_ref.chi2min);
fprintf('rw  fit: m^2 = %.2g eV^2 , chi2min = %.2g \n',FitResults.par(1),FitResults.chi2min);