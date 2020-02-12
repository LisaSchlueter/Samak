addpath(genpath('../../../Samak2.0'));
close all

run = 40257;
reffile = sprintf('Run%g.mat',run);
mtd = load(reffile);
TD      = sprintf('Run%g',run);
TimeSec = mtd.RunTime;
initrunfct = str2func(sprintf('ref_run%g',run));
Sref = initrunfct('TD',TD,'TimeSec',TimeSec,'WGTS_CD_MolPerCm2',5e17*0.9);
Sref.ComputeTBDDS;Sref.ComputeTBDIS

% Define your CM Class with your number of trials, choose SysEffect,...
myEffects = struct(...
    'BM1S','ON');
C = CovarianceMatrix('StudyObject',Sref,...
    'nTrials',0,...
    'RunFlag','OFF', 'nRuns',1,...
    'SysEffect',myEffects,...
    'BM1_RatePerSec',0.33,...
    'BM1_RedFactor',1,...
    'BM1_RefSlopeqU',16575,...
    'RecomputeFlag','ON');
C.ComputeCM_BM1S('Model','D');
C.ComputeCM;

% Save
covmatBfrac = C.MultiCovMatFrac.CM_BM1S;
save('covmatBfracRun4057.mat','covmatBfrac');

% Plot
%C.PlotCM

% Decompose
%C.DecomposeCM

% Plot
%%
figure(1)
imagesc(C.CovMatFrac);  
colormap(flipud(gray)); 
colorbar
xlabel('qU');
ylabel('qU');
title('Stacked fractional covariance matrix');
PrettyFigureFormat
publish_figure(1,'BM1S_CovMat_FT-TL3_1-Runs.mat.eps');
