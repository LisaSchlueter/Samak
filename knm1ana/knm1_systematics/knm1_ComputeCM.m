% script to demonstrate how systematic effects can be used
% KNM1 settings

%% create MultiRunAnalysis object
R = MultiRunAnalysis('RunList','KNM1',... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
    'chi2','chi2Stat',...                 % uncertainties: statistical or stat + systematic uncertainties
    'DataType','Real',...                 % can be 'Real' or 'Twin' -> Monte Carlo
    'fixPar','mNu E0 Norm Bkg',...        % free Parameter!!
    'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
    'NonPoissonScaleFactor',1.064,...     % background uncertainty are enhanced
    'minuitOpt','min ; minos',...         % technical fitting options (minuit)
    'FSDFlag','Sibille0p5eV',...               % final state distribution
    'ELossFlag','KatrinT2',...            % energy loss function
    'SysBudget',22,...                    % defines syst. uncertainties -> in GetSysErr.m;
    'DopplerEffectFlag','FSD_Knm1');      

range = 40; % eV below the endpoint
R.exclDataStart = R.GetexclDataStart(range); % set region of interest

%% default systematics
R.chi2 = 'chi2CMShape';
%R.ComputeCM; % calculate covariance matrices: no further aguments == all default sys

%% customize systematics example: only FSD systematics included
mySysEffects = struct(...
    'RF_EL','OFF',...       % Response Function(RF) EnergyLoss
    'RF_BF','OFF',...       % RF B-Fields
    'RF_RX','OFF',...       % Column Density, inel cross ection
    'FSD','ON',...          % final states
    'TASR','OFF',...        % tritium activity fluctuations
    'TCoff_RAD','OFF',...   % radiative thoretical corrections (this has to be OFF for KNM1, because included in model)
    'TCoff_OTHER','OFF',... % other theo corr.
    'DOPoff','OFF',...      % doppler effects --> OFF, because in model
    'Stack','OFF',...       % stacking / HV fluctuations
    'FPDeff','OFF');        % detector efficiency

R.ComputeCM('SysEffects',mySysEffects,...
           'BkgCM','OFF');  % background systematics (slope!) switched on/off here

%% display covariance matrices: example : FSD
R.FitCM_Obj.nTrials = 5000;
R.FitCM_Obj.ComputeCM_FSD;
R.FitCM_Obj.PlotCM;
       