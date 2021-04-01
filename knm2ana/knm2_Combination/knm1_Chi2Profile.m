% KNM-1 chi2 profile
%% settings
nFitMax  = 50;
RunList = 'KNM1';
range   = 40;         % 40eV range = 27 subruns
chi2 = 'chi2Stat';
if strcmp(chi2,'chi2Stat')
    NP = 1;
else
    NP = 1.064;
end
% Init Model Object and covariance matrix object
Real = MultiRunAnalysis('RunList',RunList,...
    'chi2',chi2,...
    'DataType','Real',...
    'fixPar','mNu E0 Norm Bkg',... free parameter
    'NonPoissonScaleFactor',NP,...
    'minuitOpt','min ; minos',...
    'FSDFlag','Sibille0p5eV',...
    'ELossFlag','KatrinT2',...
    'SysBudget',24,...
    'RadiativeFlag','ON',...
    'SynchrotronFlag','ON',...
    'AngularTFFlag','OFF',...
    'DopplerEffectFlag','OFF');
Real.exclDataStart = Real.GetexclDataStart(range);
Real.InitModelObj_Norm_BKG('RecomputeFlag','ON');
Real.ComputeCM;

%% 
Real.ResetFitResults;
Real.FitResult.par(1) = -0.8;
ScanResultsCM = Real.GetAsymFitError('Mode','Uniform','ParScanMax',1.8,'nFitMax',nFitMax,'SanityPlot','ON');

% %% plot
% Real.chi2 = 'chi2CMShape';
% pHandle = Real.PlotChi2Curve('Parameter','mNu','ScanResult',ScanResultsCM,'PlotBf','ON','HoldOn','OFF');
% %Real.chi2 = 'chi2Stat';
% %Real.PlotChi2Curve('Parameter','mNu','ScanResult',ScanResultsStat,'PlotBf','ON','HoldOn','ON');

