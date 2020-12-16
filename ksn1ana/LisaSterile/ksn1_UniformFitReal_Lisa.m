
%% create MultiRunAnalysis object
range = 40; % eV below the endpoint

R = MultiRunAnalysis('RunList','KNM1',... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
    'chi2','chi2Stat',...%CMShape',...              % uncertainties: statistical or stat + systematic uncertainties
    'DataType','Real',...                 % can be 'Real' or 'Twin' -> Monte Carlo
    'fixPar','mNu E0 Norm Bkg mnu4Sq sin2T4',...        % free Parameter!!
    'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
    'NonPoissonScaleFactor',1.064,...     % background uncertainty are enhanced
    'minuitOpt','min ; minos',...         % technical fitting options (minuit)
    'FSDFlag','SibilleFull',...           % final state distribution
    'ELossFlag','KatrinT2',...            % energy loss function
    'SysBudget',22,...
    'SynchrotronFlag','ON',...
    'AngularTFFlag','OFF',...
    'pullFlag',21,...
    'TwinBias_Q',18573.7);          
%%
R.exclDataStart = R.GetexclDataStart(range);
R.Fit
FitResult = R.FitResult;

