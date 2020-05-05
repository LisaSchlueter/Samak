% KNM1 basic analysis
% FPD - uniform

%% settings
range = 40; % eV below the endpoint
%% create MultiRunAnalysis object

R = MultiRunAnalysis('RunList','KNM1',... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
    'chi2','chi2Stat',...                 % uncertainties: statistical or stat + systematic uncertainties
    'DataType','Real',...                 % can be 'Real' or 'Twin' -> Monte Carlo
    'fixPar','mNu E0 Norm Bkg',...        % free Parameter!!
    'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
    'NonPoissonScaleFactor',1.064,...     % background uncertainty are enhanced
    'minuitOpt','min ; minos',...         % technical fitting options (minuit)
    'FSDFlag','SibilleFull',...               % final state distribution
    'ELossFlag','KatrinT2',...            % energy loss function
    'SysBudget',22,...
    'DopplerEffectFlag','FSD_Knm1',...
    'Twin_SameCDFlag','ON',...
    'Twin_SameIsotopFlag','ON');         % defines syst. uncertainties -> in GetSysErr.m;

% R.ModelObj.mnu4Sq_i=200;
R.exclDataStart = R.GetexclDataStart(range); % set region of interest

%% look at some properties
%R.RunData  % this is the stacked trititum data
%R.InitModelObj_Norm_BKG('Recompute','ON');
%R.ModelObj % this is the TBD object   --> model
% R.ModelObj.ComputeTBDDS; % the differential spectrum is calculated
% R.ModelObj.ComputeTBDIS; % the integral spectrum is calculated

%% Fit
R.Fit;

R.PrintDataStatistics
%% Display fit and save plot
% R.PlotFit('LabelFlag','FinalKNM1',...
%     'saveplot','pdf',...
%     'ErrorBarScaling',1,...
%     'YLimRes',[-2.2,2.9],...
%     'Colors','RGB',...
%     'DisplayStyle','PRL',...
%     'FitResultsFlag','OFF',...
%     'qUDisp','Rel');
