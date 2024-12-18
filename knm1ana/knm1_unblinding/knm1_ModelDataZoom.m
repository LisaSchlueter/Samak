% KNM1 Final fit Results
% Uniform Fit
% Golden Run List
% Golden Pixel List
%% settings
RunList = 'KNM1';
range   = 40;         % 40eV range = 27 subruns
% Init Model Object and covariance matrix object
Real = MultiRunAnalysis('RunList',RunList,...
    'chi2','chi2Stat',...
    'DataType','Real',...
    'fixPar','mNu E0 Norm Bkg',... free parameter
    'RadiativeFlag','ON',...
    'NonPoissonScaleFactor',1.064,...
    'minuitOpt','min ; minos',...
    'FSDFlag','Sibille0p5eV',...
    'ELossFlag','KatrinT2',...
    'SysBudget',22,...new
    'AngularTFFlag','OFF');
Real.exclDataStart = Real.GetexclDataStart(range); 
%% Get Covariance Matrix
Real.chi2='chi2CMShape';%'; % shape-only mode
Real.ComputeCM('BkgMode','SlopeFit','nTrials',1000);%('SysEffect',struct('FSD','ON','TASR','ON','Stacking','ON'),'BkgCM','ON');
Real.Fit;
%%
close all
Real.PlotDataModel_KNM1('TickDir','Out');
