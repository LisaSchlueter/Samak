% Plot Covariance Matrix for 1 single systematic effect
% systematics setting

myEffect      = 'LongPlasma';
RecomputeFlag = 'OFF';
nTrials       =  1000;
SysBudget     = 34; % 34= knm2 preliminary input

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
if strcmp(myEffect,'Bkg')
    CmArg = {'BkgCM','ON',...%  'SysEffects',struct(myEffect,'ON'),...
        'RecomputeFlag',RecomputeFlag,...
        'nTrials',nTrials,...
        'SysEffects',struct('FSD','OFF'),...
        'PlotSaveCM','ON'};
else
    CmArg = {'BkgCM','OFF',...%  'SysEffects',struct(myEffect,'ON'),...
        'RecomputeFlag',RecomputeFlag,...
        'nTrials',nTrials,...
        'SysEffects',struct(myEffect,'ON'),...
        };
end
%% init model
M = MultiRunAnalysis(RunArg{:});
M.chi2 = 'chi2CMShape';

%% label
plotdir = [getenv('SamakPath'),'knm2ana/knm2_systematics/plots/'];

if ~strcmp(myEffect,'all')
    %% calculate and plot specific covariance matrix
    M.ComputeCM(CmArg{:});
    % display and save to plots
    M.FitCM_Obj.PlotCM('qUWindowIndexMax',10,'saveplot','ON','Convergence',...
        'ON','Mode','Shape','PlotEffect',myEffect,'savedir',plotdir);
    M.FitCM_Obj.PlotCorr('qUWindowIndexMax',10,'saveplot','ON',...
        'savedir',plotdir,'savename',myEffect);
    
elseif strcmp(myEffect,'all')
    %% plot total stat + syst CM
    M.ComputeCM;
    
    % cov mat frac shape
    M.FitCM_Obj.PlotCM('qUWindowIndexMax',10,'saveplot',...
        'ON','Convergence','OFF','CovMatInput',M.FitCMFracShape,'PlotEffect','total',...
        'savedir',plotdir,'savename','FS_Uniform');
    % correlation
    M.FitCM_Obj.PlotCorr('qUWindowIndexMax',10,'saveplot',...
        'ON','CovMatInput',M.FitCMFracShape,...
        'savedir',plotdir,'savename','KNM1_FS_Uniform');
end

