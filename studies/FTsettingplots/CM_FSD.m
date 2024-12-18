addpath(genpath('../../../Samak2.0'));
close all

% Build FSD Covariance Matrix with a given TD
A = ref_FTMTD();

% Define your CM Class with your number of trials, choose SysEffect,...
myEffects = struct(...
    'TCoff','OFF',...
    'FSD','ON');
C = CovarianceMatrix('StudyObject',A,...
    'nTrials',100,...
    'RunFlag','OFF', 'nRuns',1,...
    'SysEffect',myEffects,...
    'FSDNorm_RelErr',5e-2,...
    'FSDShape_RelErr',5e-2,...
    'RecomputeFlag','ON');
C.ComputeCM_FSD;
%C.ComputeCM;

% Plot
C.PlotCM

% Decompose
C.DecomposeCM
