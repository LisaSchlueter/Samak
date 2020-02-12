% Script to Develop, Test and Display 
% Covariance matrices for multi-ring fit
% taking ring to ring correlation into account
% Lisa, Oct 2019

%% Set up multi ring model
RunList = 'KNM2_Prompt';
exclDataStart= 14;
chi2 = 'chi2Stat';
RingMerge = 'Half'; % 'Full == 4 rings

% read data and set up model
CommonArg = {'RunList',RunList,...
    'chi2','chi2Stat','DataType','Real',...
    'exclDataStart',exclDataStart,...
    'fixPar','mNu Norm E0 Bkg',...
    'RadiativeFlag','ON',...
    'minuitOpt','min ; minos',...
    'FSDFlag','Sibille0p5eV',...
    'ELossFlag','KatrinT2',...
    'SysBudget',22,...
    'AnaFlag','Ring',...
    'RingMerge',RingMerge,...
    'chi2',chi2};

R = MultiRunAnalysis(CommonArg{:});
%R.Fit;

%% Compute multi ring covariance matrix
%R.ComputeCM('SysEffects',struct('FSD','ON'),'BkgCM','OFF','nTrials',1000);
%R.ComputeCM('SysEffects',struct('TC','ON'),'BkgCM','OFF');
%R.ComputeCM('SysEffects',struct('TC','OFF'),'BkgCM','ON');
%R.ComputeCM('SysEffects',struct('Stack','ON'),'BkgCM','OFF');
%R.ComputeCM('SysEffects',struct('TASR','ON'),'BkgCM','OFF');

%% 
%R.chi2='chi2CMShape';
%R.Fit