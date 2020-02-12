addpath(genpath('../../../Samak2.0'));
close all

% 1. Define your Study Obect (TBD Object)
A = ref_TASR();

% 2. Define your CM Class with your number of trials, choose SysEffect,...
myEffects = struct(...
    'RF_EL','ON',... % Response Function(RF) EnergyLoss
    'RF_BF','ON',... % RF B-Fields
    'RF_RX','ON',... % RF Column Density, Cross Section
    'TASR','ON',... % Tritium Column Density Stack Runs);
    'FSD','ON');

C = CovarianceMatrix('StudyObject',A, 'nTrials',1000,...
    'WGTS_TASR_RelErr',0.01,...
    'RunFlag','OFF', 'RecomputeFlag','ON');

% 3. Compute TASR Covariance Matrix
C.RecomputeFlag = 'ON'; %Enforces Recomputing the Matrix
C.ComputeCM_TASR; C.ComputeCM_RF;  C.ComputeFracCM;  


%% Plot Stat
StatFrac = bsxfun(@rdivide,diag(A.TBDIS),A.TBDIS);  %divides 1st row of CovMat by TBDIS(1), etc...
StatFrac = bsxfun(@rdivide,StatFrac,A.TBDIS'); %divides columnwise
                    
figure(1)
imagesc(StatFrac);  
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

