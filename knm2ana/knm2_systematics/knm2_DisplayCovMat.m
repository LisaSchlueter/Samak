% Plot Covariance Matrix for 1 single systematic effect
% systematics setting

myEffect      = 'LongPlasma';
RecomputeFlag = 'OFF';
nTrials       =  1000; 
SysBudget     = 33; % 33= knm2 preliminary input

% model setting
RunList   = 'KNM2_Prompt';
AnaFlag   = 'StackPixel'; % FPD segmentation

E0 = knm2FS_GetE0Twins('SanityPlot','OFF');

RunArg = {'FSDFlag','BlindingKNM2',...
    'ELossFlag','KatrinT2',...
    'AnaFlag',AnaFlag,...
    'chi2','chi2Stat',...      % switch on later by hand
    'RunList',RunList,...
    'fixPar','E0 Norm Bkg',... % free parameter
    'DataType','Twin',...      % MC twins
    'TwinBias_Q',E0,...  % twin endpoint
    'SysBudget',SysBudget,...
    'RingMerge','Full',...
    'NonPoissonScaleFactor',1.1121};

CmArg = {'BkgCM','OFF',...%  'SysEffects',struct(myEffect,'ON'),...
    'RecomputeFlag',RecomputeFlag,...
    'nTrials',nTrials,...
    'SysEffects',struct(myEffect,'ON'),...
    };
%% init model
M = MultiRunAnalysis(RunArg{:});
M.chi2 = 'chi2CMShape';
M.ComputeCM;
%% calculate covariance matrix
M.ComputeCM(CmArg{:});
%% display and save to plots
M.FitCM_Obj.PlotCM('qUWindowIndexMax',10,'saveplot','ON','Convergence','ON');





