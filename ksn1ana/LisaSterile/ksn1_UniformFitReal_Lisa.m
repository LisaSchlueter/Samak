
range = 95; % eV below the endpoint
%% create MultiRunAnalysis object

R = MultiRunAnalysis('RunList','KNM1',... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
    'chi2','chi2Stat',...%CMShape',...              % uncertainties: statistical or stat + systematic uncertainties
    'DataType','Twin',...                 % can be 'Real' or 'Twin' -> Monte Carlo
    'fixPar','E0 Norm Bkg',...        % free Parameter!!
    'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
    'NonPoissonScaleFactor',1.064,...     % background uncertainty are enhanced
    'minuitOpt','min ; minos',...         % technical fitting options (minuit)
    'FSDFlag','SibilleFull',...           % final state distribution
    'ELossFlag','KatrinT2',...            % energy loss function
    'SysBudget',22,...
    'SynchrotronFlag','ON',...
    'AngularTFFlag','OFF',...
    'pullFlag',9,...
    'TwinBias_Q',18573.7);          
%%
R.exclDataStart = R.GetexclDataStart(95);
R.Fit
FitResults_95 = R.FitResult;
%%
R.exclDataStart =12;% used for twinsR.GetexclDataStart(40);
R.Fit
FitResults_12 = R.FitResult;
%%
% R.ModelObj.mnu4Sq_i=200;
R.exclDataStart = R.GetexclDataStart(range); % set region of interest
m4_bf     = 56.234; % eV
sin2t4_bf = 0.0133;
R.ModelObj.SetFitBiasSterile(m4_bf ^2,sin2t4_bf); % fix sterile parameter to minimal chi2 from grid search
R.Fit;
FitResult_Sterile = R.FitResult;
%%

R.ModelObj.SetFitBiasSterile(0,0); % fix sterile parameter to minimal chi2 from grid search
R.Fit;
FitResult_Null = R.FitResult;
%%

fprintf('sterile (best-fit) : chi2 = %.2f (%.0f dof), p-values = %.2f \n',...
    FitResult_Sterile.chi2min,FitResult_Sterile.dof,1-chi2cdf(FitResult_Sterile.chi2min,FitResult_Sterile.dof));
fprintf('no sterile         : chi2 = %.2f (%.0f dof), p-values = %.2f \n',...
    FitResult_Null.chi2min,FitResult_Null.dof,1-chi2cdf(FitResult_Null.chi2min,FitResult_Null.dof));



%%
