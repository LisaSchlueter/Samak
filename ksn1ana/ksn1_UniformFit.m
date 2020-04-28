% KNM1 basic analysis
% FPD - uniform

%% settings
range = 95; % eV below the endpoint
%% create MultiRunAnalysis object

R = MultiRunAnalysis('RunList','KNM1',... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
    'chi2','chi2CMShape',...                 % uncertainties: statistical or stat + systematic uncertainties
    'DataType','Real',...                 % can be 'Real' or 'Twin' -> Monte Carlo
    'fixPar','mNu E0 Norm Bkg mnu4Sq sin2T4',...        % free Parameter!!
    'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
    'NonPoissonScaleFactor',1.064,...     % background uncertainty are enhanced
    'minuitOpt','min ; minos',...         % technical fitting options (minuit)
    'FSDFlag','SibilleFull',...           % final state distribution
    'ELossFlag','KatrinT2',...            % energy loss function
    'SysBudget',22,...
    'DopplerEffectFlag','FSD_Knm1',...
    'Twin_SameCDFlag','ON',...
    'Twin_SameIsotopFlag','ON',...
    'SynchrotronFlag','OFF',...
    'AngularTFFlag','OFF',...
    'pullFlag',9);          

% R.ModelObj.mnu4Sq_i=200;
R.exclDataStart = R.GetexclDataStart(range); % set region of interest

%% Fit
R.Fit;

%% Display fit and save plot
R.PlotFit('LabelFlag','data',...
    'saveplot','pdf',...
    'ErrorBarScaling',1,...
    'YLimRes',[-2.2,2.9],...
    'Colors','RGB',...
    'DisplayStyle','PRL',...
    'FitResultsFlag','OFF',...
    'qUDisp','Rel');


% ===============================================
%   m^2       = -0.985649 +/- 0.48558 eV^2
%   m         = 0 +/- 0.696836 eV
%  - - - - - - - - - - - - - - - - - - - - - - - 
%   mnu4Sq    = 0.0196248 +/- 0.376255 eV
%   sin2T4    = 0.0107855 +/- 0.0418495 eV
%  - - - - - - - - - - - - - - - - - - - - - - - 
%   (E0)eff   = 18573.7164 +/- 0.0206743 eV
%  - - - - - - - - - - - - - - - - - - - - - - - 
%   B         = 292.336 +/- 0.666961 mcps (117 pixels) 
%  - - - - - - - - - - - - - - - - - - - - - - - 
%   N         = 0.835759 +/- 0.00031119
%  - - - - - - - - - - - - - - - - - - - - - - - 
%   chi2/dof  = 41.2665/33
%   p-value   = 0.152975
% ===============================================
% ===============================================
%   m^2       = -0.554448 +/- 0.647021 eV^2
%   m         = 0 +/- 0.804376 eV
%  - - - - - - - - - - - - - - - - - - - - - - - 
%   mnu4Sq    = 0.019237 +/- 0.343557 eV
%   sin2T4    = 0.0131271 +/- 0.0751944 eV
%  - - - - - - - - - - - - - - - - - - - - - - - 
%   (E0)eff   = 18573.7315 +/- 0.0695626 eV
%  - - - - - - - - - - - - - - - - - - - - - - - 
%   B         = 292.458 +/- 0.71799 mcps (117 pixels) 
%  - - - - - - - - - - - - - - - - - - - - - - - 
%   N         = 0.835707 +/- 0.00131674
%  - - - - - - - - - - - - - - - - - - - - - - - 
%   chi2/dof  = 32.8646/33
%   p-value   = 0.473882
% ===============================================
