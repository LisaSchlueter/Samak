addpath(genpath('../../../Samak2.0'));
close all

% 1. Define your Study Obect (TBD Object)
A  = ref_FTMTD();
AD = ref_FTMTD('DopplerEffectFlag','matConv');

% 2. Define your CM Class with your number of trials, choose SysEffect,...
myEffects = struct(...
    'STAT','ON',...
    'RF_EL','ON',...  % Response Function(RF) EnergyLoss
    'RF_BF','ON',...  % RF B-Fields
    'RF_RX','ON',...  % RF Column Density, Cross Section
    'TASR','ON',...   % Tritium Activity Stack Runs);
    'FSD','ON',...     % FSD Noralization
    'DOPoff','OFF',...  % Doppler (if switched off)
    'TCoff_RAD','ON',...    % theretical corrections (if switched off)
    'TCoff_OTHER','ON',... % theretical corrections (if switched off)
    'BM1S','OFF');      % Background

C = CovarianceMatrix('StudyObject',A, 'AltObject', AD,...
    'nTrials',1000,...
    'RunFlag','OFF', 'nRuns',1,...
    'SysEffect',myEffects,...
    'RecomputeFlag','OFF',...
    'MACE_Ba_T_RelErr',1e-2,...
    'MACE_Bmax_T_RelErr',1e-2,...
    'WGTS_B_T_RelErr',1e-2,...
    'WGTS_CD_MolPerCm2_RelErr',0.1,...
    'ISXsection_RelErr',1e-2,...
    'WGTS_TASR_RelErr',1e-2,...
    'FSDNorm_RelErr',3e-2,...
    'FSDShape_RelErr',3e-2,...
    'BM1_RatePerSec',0.3,...
    'BM1_RedFactor',0.5,...
    'BM1_RefSlopeqU',16375);

C.ComputeCM;  
C.PlotStack('plotStat','ON');

%%
figure(1)
imagesc(C.CovMatFrac);  
colormap(flipud(gray)); 
colorbar
xlabel('qU');
ylabel('qU');
title('Stacked fractional covariance matrix');
PrettyFigureFormat
publish_figure(1,'CovMat_FT-TL4_1000-Trials_-RF_EL-RF_BF-RF_RX-TASR-FSD-TCoff_RAD-TCoff_OTHER.eps');

return;

%% Plot Stat

figure(1)
imagesc(C.MultiCovMatFrac.CM_STAT);  
colormap(flipud(gray)); 
colorbar
xlabel('qU');
ylabel('qU');
title('Statistical fluctuations');
PrettyFigureFormat
publish_figure(1,'stat.eps');


%% Plot CM_TASR
figure(2)
imagesc(C.MultiCovMatFrac.CM_TASR);  
colormap(flipud(gray)); 
colorbar
xlabel('qU');
ylabel('qU');
title('4% subrun activity fluctuations');
PrettyFigureFormat
publish_figure(2,'tasr.eps');


%% Plot RF
figure(3)
imagesc(C.MultiCovMatFrac.CM_RF);  
colormap(flipud(gray)); 
colorbar
xlabel('qU');
ylabel('qU');
title('Response function');
PrettyFigureFormat
publish_figure(3,'rf.eps');


%% Plot RF+Stat
figure(4)
imagesc(C.MultiCovMatFrac.CM_RF+StatFrac);  
colormap(flipud(gray)); 
colorbar
xlabel('qU');
ylabel('qU');
title('Stat + Response function');
PrettyFigureFormat
publish_figure(4,'statrf.eps');

%% Plot TASR+Stat
figure(5)
imagesc(C.MultiCovMatFrac.CM_TASR+StatFrac);  
colormap(flipud(gray)); 
colorbar
xlabel('qU');
ylabel('qU');
title('Stat + Activity Fluctuation');
PrettyFigureFormat
publish_figure(5,'stattasr.eps');

%% Plot TASR+Stat+RF
figure(6)
imagesc(C.MultiCovMatFrac.CM_TASR+C.MultiCovMatFrac.CM_RF+StatFrac);  
colormap(flipud(gray)); 
colorbar
xlabel('qU');
ylabel('qU');
title('Stat + Response function + Activity Fluctuation');
PrettyFigureFormat
publish_figure(6,'stattasrrf.eps');

C.DecomposeCM;  % Plot Decomposed CM

