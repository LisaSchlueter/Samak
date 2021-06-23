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
    'FSDFlag','SibilleFull',...
    'ELossFlag','KatrinT2',...
    'SysBudget',22,...
    'AngularTFFlag','OFF');
Real.exclDataStart = Real.GetexclDataStart(range); 
%% Get Covariance Matrix
Real.chi2='chi2CMShape';%'; % shape-only mode
Real.ComputeCM('BkgMode','SlopeFit','nTrials',1000);%('SysEffect',struct('FSD','ON','TASR','ON','Stacking','ON'),'BkgCM','ON');
Real.ModelObj.mnuSq_i = -0.98;
Real.Fit;
%% Fit
%Real.Fit('CATS','OFF');
FitCMShape= Real.FitCMShape;
%% Display fit and save plot
Real.FitCMShape = FitCMShape.*1.08;
Real.PlotFit('LabelFlag','FinalKNM1',...
    'saveplot','pdf',...
    'ErrorBarScaling',50,...
    'YLimRes',[-2.2,2.9],...
    'Colors','RGB',...
    'DisplayStyle','PRL',...
    'FitResultsFlag','OFF',...
    'qUDisp','Abs',...
    'TickDir','Out');
