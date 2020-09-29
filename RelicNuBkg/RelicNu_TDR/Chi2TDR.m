%% settings
range = 30; % eV below the endpoint
Netabins = 10;
etarange = 9;
etafactor = 3; %max(eta) = etafactor*10^(etarange+1)
initfile=@ref_RelicNuBkg_DesignReport;

R = RunAnalysis('RunNr',1,... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
    'FakeInitFile',initfile,...
    'chi2','chi2Stat',...                 % uncertainties: statistical or stat + systematic uncertainties
    'DataType','Fake',...                 % can be 'Real' or 'Twin' -> Monte Carlo
    'fixPar','mNu E0 Norm Bkg',...        % free Parameter!!
    'NonPoissonScaleFactor',1,...         % background uncertainty are enhanced
    'minuitOpt','min ; minos',...         % technical fitting options (minuit)
    'FSDFlag','Sibille0p5eV',...          % final state distribution                        !!check ob initfile hier überschrieben wird
    'ELossFlag','KatrinT2',...            % energy loss function
    'SysBudget',22,...                    % defines syst. uncertainties -> in GetSysErr.m;
    'DopplerEffectFlag','FSD',...
    'SynchrotronFlag','OFF',...
    'AngularTFFlag','OFF');

R.exclDataStart = 1; % set region of interest

%R.InitModelObj_Norm_BKG('Recompute','ON');

Chi2 = 1:Netabins;
mnu  = 1:Netabins;
E0   = 1:Netabins;
Bkg  = 1:Netabins;

for i=1:Netabins
   R.ModelObj.eta = (i-1)*((etafactor*10^(etarange))/(Netabins-1));
   R.ModelObj.ComputeNormFactorTBDDS;
   R.ModelObj.ComputeTBDDS;
   R.ModelObj.ComputeTBDIS;
   %% Fit
   R.Fit;
   %R.PlotFit;
   Chi2(i)=R.FitResult.chi2min;
   mnu(i)=R.ModelObj.mnuSq_i+R.FitResult.par(1);
   E0(i)=R.ModelObj.Q_i+R.FitResult.par(2);
   Bkg(i)=R.ModelObj.BKG_RateSec_i+R.FitResult.par(3);
end

save('RelicChi2Scan_TDR1.mat','Chi2','Netabins','etafactor','etarange','mnu','E0','Bkg');