% KNM1 basic analysis
% FPD - uniform

%% create MultiRunAnalysis object
R = MultiRunAnalysis('RunList','KNM1',...        % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
    'chi2','chi2CMShape',...                     % uncertainties: statistical or stat + systematic uncertainties
    'DataType','Real',...                        % can be 'Real' or 'Twin' -> Monte Carlo
    'fixPar','mNu E0 Norm Bkg mnu4Sq sin2T4',...                   % free Parameter!!
    'RadiativeFlag','ON',...                     % theoretical radiative corrections applied in model
    'NonPoissonScaleFactor',1.064,...            % background uncertainty are enhanced
    'minuitOpt','min ; minos',...               % technical fitting options (minuit)
    'FSDFlag','SibilleFull',...                  % final state distribution
    'ELossFlag','KatrinT2',...                   % energy loss function
    'SysBudget',22);                             % defines syst. uncertainties -> in GetSysErr.m;

R.pullFlag=1;
%R.ModelObj.mnuSq_i=-0.956847;
% R.i_Q=0;
% R.i_B=[];
% R.i_N=[];
% R.ModelObj.mnu4Sq_i=0;
% R.ModelObj.sin2T4_i=0;
R.exclDataStart = R.GetexclDataStart(90); % set region of interest
% R.ModelObj.PlotFSD
%% look at some properties
R.RunData  % this is the stacked trititum data
%R.ModelObj % this is the TBD object   --> model
% R.ModelObj.ComputeTBDDS; % the differential spectrum is calculated
% R.ModelObj.ComputeTBDIS; % the integral spectrum is calculated

%% Fit
R.Fit;
X = R.FitResult.chi2min-0.0101
%% Display fit and save plot
% R.PlotFit('LabelFlag','FinalKNM1',...
%     'saveplot','pdf',...
%     'ErrorBarScaling',1,...
%     'YLimRes',[-2.2,2.9],...
%     'Colors','RGB',...
%     'DisplayStyle','PRL',...
%     'FitResultsFlag','OFF',...
%     'qUDisp','Rel');