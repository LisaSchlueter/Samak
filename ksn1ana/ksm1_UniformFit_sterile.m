% KNM1 basic analysis
% FPD - uniform

%% settings
range = 90; % eV below the endpoint
%% create MultiRunAnalysis object
R = MultiRunAnalysis('RunList','KNM1',...        % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
    'chi2','chi2CMShape',...                     % uncertainties: statistical or stat + systematic uncertainties
    'DataType','Twin',...                        % can be 'Real' or 'Twin' -> Monte Carlo
    'fixPar','E0 Norm Bkg',...                   % free Parameter!!
    'RadiativeFlag','ON',...                     % theoretical radiative corrections applied in model
    'NonPoissonScaleFactor',1.064,...            % background uncertainty are enhanced
    'minuitOpt','min ; migrad',...               % technical fitting options (minuit)
    'FSDFlag','Sibille0p5eV',...                 % final state distribution
    'ELossFlag','KatrinT2',...                   % energy loss function
    'SysBudget',22);                             % defines syst. uncertainties -> in GetSysErr.m;

%R.pullFlag=9;
%R.ModelObj.mnuSq_i=-0.956847;
% R.ModelObj.mnu4Sq_i=200;
% R.ModelObj.sin2T4_i=0.07;
R.exclDataStart = R.GetexclDataStart(40); % set region of interest

%% look at some properties
R.RunData  % this is the stacked trititum data
%R.ModelObj % this is the TBD object   --> model
% R.ModelObj.ComputeTBDDS; % the differential spectrum is calculated
% R.ModelObj.ComputeTBDIS; % the integral spectrum is calculated

%% Fit
R.Fit;

%% Display fit and save plot
R.PlotFit('LabelFlag','FinalKNM1',...
    'saveplot','pdf',...
    'ErrorBarScaling',1,...
    'YLimRes',[-2.2,2.9],...
    'Colors','RGB',...
    'DisplayStyle','PRL',...
    'FitResultsFlag','OFF',...
    'qUDisp','Rel');