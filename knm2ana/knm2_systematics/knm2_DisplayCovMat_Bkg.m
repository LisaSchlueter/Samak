% Plot Covariance Matrix for background (qU) sloep systematic effect
% systematics setting
RecomputeFlag = 'OFF';
nTrials       = 100; 
SysBudget     = 31; % 31= knm2 preliminary input

% model setting
RunList   = 'KNM2_RW1';
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

CmArg = {'BkgCM','ON',...
    'SysEffects',struct('BkgShape','ON'),...
    'RecomputeFlag',RecomputeFlag,...
    'nTrials',nTrials,...
    };
%% init model
M = MultiRunAnalysis(RunArg{:});

%% calculate covariance matrix
M.chi2 = 'chi2CMShape';
M.ComputeCM(CmArg{:});

%% display and save to plots
M.FitCM_Obj.PlotCM('qUWindowIndexMax',10,'saveplot','ON');

%% 
M.InitModelObj_Norm_BKG;
M.FitCM_Obj.ComputeCM_Background('Display','ON');




