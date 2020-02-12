% script to do estimate upper DT uncertainty
%for First Tritium
% Plan: 1. simulate spectrum, where with different DT fluctuations
%       2. add stat flutuations
%       % fit with statistics: how good is chi2?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
% General settings
myRunList = {'StackCD100all','StackCD100up','StackCD100down'}; 
ringCutFlag = 'ex2';
%Fit settings
fixPar = '1 5 6';
chi2Flag = 'chi2Stat';
% Simulation settings
StatFluct = 'OFF'; % add stat fluctuation 
SysShift = 'ON';  % add systematic shift
nuMass = 0;         %Neutrino mass of simulation model
nSim = 100;         %Number of simulated spectra
mySysEffectsSim =  struct('TASR','ON'); % sys effects in MC data
WGTS_TASR_RelErr = [0.001 0.003 0.005 0.01 0.02 0.03 0.04];
nTrials = 1000;
RecomputeFlag = 'OFF';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize RunAnalysis Model
A = MultiRunAnalysis('RunList',myRunList{1},'AnaFlag','StackPixel',...
    'fixPar',fixPar,'pullFlag',3);

% Init other variables
TBDISsim       = zeros(A.ModelObj.nqU,nSim,numel(WGTS_TASR_RelErr));
FitResults_par = zeros(6,nSim,numel(WGTS_TASR_RelErr));
FitResults_err = zeros(6,nSim,numel(WGTS_TASR_RelErr));
Fitchi2min     = zeros(nSim,numel(WGTS_TASR_RelErr));
meanChi2min    = zeros(numel(WGTS_TASR_RelErr),1);
stdChi2min     = zeros(numel(WGTS_TASR_RelErr),1);

% Initialize Simulation 
S = ref_RunAnalysis(myRunList{1},ringCutFlag,'','ISCS','Theory',...
                    'recomputeRF','OFF','mnuSq_i',nuMass);              
S.ComputeTBDDS;
S.ComputeTBDIS;
TBDIS_i = S.TBDIS;

% Add Fluctuation: Stat and/or Sys
for ii=1:numel(WGTS_TASR_RelErr)
if strcmp(SysShift,'ON')
    A.ComputeCM('SysEffects',mySysEffectsSim,'nTrials',nTrials,...
        'RecomputeFlag',RecomputeFlag,'InitNormFit','OFF',...
        'StackCM','OFF','WGTS_TASR_RelErr',WGTS_TASR_RelErr(ii));
    
    CMFrac = A.FitCM_Obj.MultiCovMatFrac.CM_TASR;
    TBDIS_i = mvnrnd(S.TBDIS,S.TBDIS.*CMFrac.*S.TBDIS',nSim)';
end
if strcmp(StatFluct,'ON')
    TBDISsim(:,:,ii) = TBDIS_i + repmat(sqrt(S.TBDIS),1,nSim).*randn(S.nqU,nSim);
elseif strcmp(StatFluct,'OFF') && strcmp(SysShift,'ON')
    TBDISsim(:,:,ii) = TBDIS_i;
elseif strcmp(StatFluct,'OFF') && strcmp(SysShift,'OFF')
    TBDISsim(:,:,ii) = repmat(TBDIS_i,1,nSim);
end

% Get Covariance Matrix for Fit (if needed)
A.chi2 = chi2Flag;
if ~strcmp(A.chi2,'chi2Stat')
    % check if mySysEffects and mySysEffectsSim  are the same
    names = fieldnames(mySysEffects);
    tmp =0;
    for i=1:numel(names) 
        if  strcmp(mySysEffects.(names{i}),mySysEffectsSim.(names{i}))
            tmp = tmp +1;
        end
    end
    % if not the same
    if tmp~=numel(names) || strcmp('SysShift','OFF')
        A.ComputeCM('SysEffects',mySysEffects,'nTrials',nTrials,...
            'RecomputeFlag',RecomputeFlag,'InitNormFit','OFF',...
            'StackCM',StackCM);
    end
end

% Do the Fit
A.RunData.qU    = S.qU;
for i=1:nSim
A.RunData.TBDIS = TBDISsim(:,i,ii);
A.Fit;
FitResults_par(:,i,ii) = A.FitResult.par;
FitResults_err(:,i,ii) = A.FitResult.err;
Fitchi2min(i,ii)       = A.FitResult.chi2min;
end
meanChi2min(ii) = mean(Fitchi2min(:,ii));
stdChi2min(ii)  = std(Fitchi2min(:,ii));
end

%% Plot
errorbar(WGTS_TASR_RelErr*100,meanChi2min,stdChi2min,'--o','LineWidth',3);
xlabel('Tritium Activity Fluctuation Uncertainty (%)');
ylabel('\chi²');
PrettyFigureFormat;
