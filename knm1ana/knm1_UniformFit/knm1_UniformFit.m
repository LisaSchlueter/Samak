% KNM1 basic analysis
% FPD - uniform

%% settings
range = 40; % eV below the endpoint
%% create MultiRunAnalysis object

R = MultiRunAnalysis('RunList','KNM1',... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
    'chi2','chi2CMShape',...                 % uncertainties: statistical or stat + systematic uncertainties
    'DataType','Real',...                 % can be 'Real' or 'Twin' -> Monte Carlo
    'fixPar','mNu E0 Norm Bkg',...        % free Parameter!!
    'NonPoissonScaleFactor',1.064,...     % background uncertainty are enhanced
    'minuitOpt','min ; minos',...         % technical fitting options (minuit)
    'FSDFlag','SibilleFull',...          % final state distribution
    'ELossFlag','KatrinT2',...            % energy loss function
    'SysBudget',22,...                    % defines syst. uncertainties -> in GetSysErr.m;
    'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
    'DopplerEffectFlag','OFF',...
    'SynchrotronFlag','ON',...
    'AngularTFFlag','OFF');

%%
R.exclDataStart = R.GetexclDataStart(range); % set region of interest
R.InitModelObj_Norm_BKG('RecomputeFlag','ON');
R.ComputeCM;
%% Fit
R.Fit;
%% cov mats tmp
savedir = [getenv('SamakPath'),'knm1ana/knm1_systematics/knm1_covmats/'];%inputs/CovMat/Background/CM/'];
filenameTot = sprintf('%sknm1_TotalFitCM.mat',savedir);
qU = R.ModelObj.qU;
Signal = R.ModelObj.TBDIS - R.ModelObj.TimeSec .* R.ModelObj.qUfrac .* R.ModelObj.BKG_RateSec;
Signal(R.RunData.qU>R.ModelObj.Q) = 0;
Background      = R.ModelObj.TimeSec .* R.ModelObj.qUfrac .* R.ModelObj.BKG_RateSec;

FitCMFracShape = R.FitCMFracShape ;
FitCMShape     = R.FitCMShape ;
FitCM          = R.FitCM;
FitCMFrac      = R.FitCMFrac;

save(filenameTot,'qU','Signal','Background','FitCM','FitCMFrac','FitCMShape','FitCMFracShape');
%% cov mats tmp
savedir = [getenv('SamakPath'),'knm1ana/knm1_systematics/knm1_covmats/'];%inputs/CovMat/Background/CM/'];
filenameStat = sprintf('%sknm1_NonPoisCM.mat',savedir);

% Seperate Signal and Background
qU = R.ModelObj.qU;
Signal = R.ModelObj.TBDIS - R.ModelObj.TimeSec .* R.ModelObj.qUfrac .* R.ModelObj.BKG_RateSec;
Signal(R.RunData.qU>R.ModelObj.Q) = 0;
SignalE         = sqrt(Signal);
%            backgroundE     = R.NonPoissonScaleFactor * sqrt(R.ModelObj.TimeSec .* R.ModelObj.qUfrac .* R.ModelObj.BKG_RateSec);%R.ModelObj.BKG_RateSec);
Background      = R.ModelObj.TimeSec .* R.ModelObj.qUfrac .* R.ModelObj.BKG_RateSec;
BackgroundE     = sqrt(R.ModelObj.TimeSec .* R.ModelObj.qUfrac .* R.ModelObj.BKG_RateSec .*R.NonPoissonScaleFactor.^2);%R.ModelObj.BKG_RateSec);

% Build statsitical covariance matrix
NP_CM       = diag(BackgroundE.^2);
Stat_CM     = diag(SignalE.^2);
StatNP_CM   = diag((SignalE.^2 + BackgroundE.^2));

StatNP_CMFrac  = diag(1./(SignalE.^2 + BackgroundE.^2));
NP_CMFrac  = diag(1./(BackgroundE.^2));
Stat_CMFrac  = diag(1./(SignalE.^2));

save(filenameStat,...
    'qU',...
    'Signal','SignalE','Background','BackgroundE',...
    'NP_CM','Stat_CM','StatNP_CM',...
    'NP_CMFrac','Stat_CMFrac','StatNP_CMFrac');
%% look at some properties
%R.RunData  % this is the stacked trititum data
%R.InitModelObj_Norm_BKG('Recompute','ON');
%R.ModelObj % this is the TBD object   --> model
% R.ModelObj.ComputeTBDDS; % the differential spectrum is calculated
% R.ModelObj.ComputeTBDIS; % the integral spectrum is calculated

return
%% Display fit and save plot
R.PlotFit('LabelFlag','data',...
    'saveplot','pdf',...
    'ErrorBarScaling',1,...
    'YLimRes',[-2.2,2.9],...
    'Colors','RGB',...
    'DisplayStyle','PRL',...
    'FitResultsFlag','OFF',...
    'qUDisp','Rel');
