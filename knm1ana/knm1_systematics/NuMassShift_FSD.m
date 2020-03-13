% KNM1 neutrino mass bias from stacking

%% settings
range = 40; % eV below the endpoint
%% create MultiRunAnalysis object
R = MultiRunAnalysis('RunList','KNM1',... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
    'chi2','chi2Stat',...                 % uncertainties: statistical or stat + systematic uncertainties
    'DataType','Twin',...                 % can be 'Real' or 'Twin' -> Monte Carlo
    'fixPar','mNu E0 Norm Bkg',...        % free Parameter!!
    'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
    'NonPoissonScaleFactor',1.064,...     % background uncertainty are enhanced
    'minuitOpt','min ; minos',...         % technical fitting options (minuit)
    'FSDFlag','Sibille0p5eV',...          % final state distribution
    'ELossFlag','KatrinT2',...            % energy loss function
    'SysBudget',22,...
    'DopplerEffectFlag','FSD_Knm1');       % defines syst. uncertainties -> in GetSysErr.m;

R.exclDataStart = R.GetexclDataStart(range); % set region of interest

%% Fit
R.Fit;

%% result
fprintf('Stacking neutrino mass bias = %.4f eV^2 \n',R.FitResult.par(1));
fprintf('Stacking fit endpoint  bias = %.4f eV \n',R.FitResult.par(2)+R.ModelObj.Q_i-R.TwinBias_Q);

