%% settings
range = 30; % eV below the endpoint
Netabins = 11;
etarange = 10;
etafactor = 2; %max(eta) = etafactor*10^(etarange)
%% create MultiRunAnalysis object

R = MultiRunAnalysis('RunList','KNM1',... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
    'chi2','chi2Stat',...                 % uncertainties: statistical or stat + systematic uncertainties
    'DataType','Twin',...                 % can be 'Real' or 'Twin' -> Monte Carlo
    'fixPar','mnu E0 Norm Bkg',...        % free Parameter!!
    'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
    'NonPoissonScaleFactor',1.064,...     % background uncertainty are enhanced
    'minuitOpt','min ; minos',...         % technical fitting options (minuit)
    'FSDFlag','Sibille0p5eV',...          % final state distribution
    'ELossFlag','KatrinT2',...            % energy loss function
    'SysBudget',22,...                    % defines syst. uncertainties -> in GetSysErr.m;
    'DopplerEffectFlag','FSD',...
    'Twin_SameCDFlag','OFF',...
    'Twin_SameIsotopFlag','OFF',...
    'SynchrotronFlag','ON',...
    'AngularTFFlag','OFF',...
    'TwinBias_Q',18573.73);

R.exclDataStart = R.GetexclDataStart(range); % set region of interest

R.InitModelObj_Norm_BKG('Recompute','ON');

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

save('RelicChi2Scan.mat','Chi2','Netabins','etafactor','etarange','mnu','E0','Bkg');

%% Display fit and save plot
%R.PlotFit('LabelFlag','data',...
%    'saveplot','pdf',...
%    'ErrorBarScaling',1,...
%    'YLimRes',[-2.2,2.9],...
%    'Colors','RGB',...
%    'DisplayStyle','PRL',...
%    'FitResultsFlag','OFF',...
%    'qUDisp','Rel');