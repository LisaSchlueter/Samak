% Plot Covariance Matrix for 1 single systematic effect
% systematics setting

myEffect      = 'RF';
RecomputeFlag = 'OFF';
nTrials       = 1000; 
SysBudget     = 31; % 31= knm2 preliminary input

% model setting
RunList   = 'KNM2_RW2';
AnaFlag   = 'StackPixel'; % FPD segmentation

RunArg = {'FSDFlag','BlindingKNM2',...
    'ELossFlag','KatrinT2',...
    'AnaFlag',AnaFlag,...
    'chi2','chi2Stat',...      % switch on later by hand
    'RunList',RunList,...
    'fixPar','E0 Norm Bkg',... % free parameter
    'DataType','Twin',...      % MC twins
    'TwinBias_Q',18573.70,...  % twin endpoint
    'SysBudget',SysBudget,...
    'RingMerge','Full'};

CmArg = {'BkgCM','OFF',...%  'SysEffects',struct(myEffect,'ON'),...
    'RecomputeFlag',RecomputeFlag,...
    'nTrials',nTrials,...
    'SysEffects',struct(myEffect,'ON'),...
    };
%% init model
M = MultiRunAnalysis(RunArg{:});

%% calculate covariance matrix
M.chi2 = 'chi2CMShape';
M.ComputeCM(CmArg{:});

%% display and save to plots
M.FitCM_Obj.PlotCM('qUWindowIndexMax',10,'saveplot','ON');





