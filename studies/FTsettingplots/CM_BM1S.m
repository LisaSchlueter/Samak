addpath(genpath('../../../Samak2.0'));
close all

% Build FSD Covariance Matrix with a given TD
A = ref_run40257();

% Define your CM Class with your number of trials, choose SysEffect,...
myEffects = struct(...
    'BM1S','ON');
C = CovarianceMatrix('StudyObject',A,...
    'nTrials',0,...
    'RunFlag','OFF', 'nRuns',1,...
    'SysEffect',myEffects,...
    'BM1_RatePerSec',3.364635822262940779e-01,...
    'BM1_RedFactor',0.2,...
    'BM1_RefSlopeqU',16575,...
    'RecomputeFlag','ON');
C.ComputeCM_BM1S('Model','D');
%C.ComputeCM;

% Plot
C.PlotCM

% Decompose
C.DecomposeCM

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
