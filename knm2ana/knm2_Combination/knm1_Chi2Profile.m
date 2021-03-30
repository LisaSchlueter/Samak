% KNM-1 chi2 profile
%% settings
nFitMax  = 50;
RunList = 'KNM1';
range   = 40;         % 40eV range = 27 subruns
chi2 = 'chi2CMShape';
% Init Model Object and covariance matrix object
Real = MultiRunAnalysis('RunList',RunList,...
    'chi2',chi2,...
    'DataType','Real',...
    'fixPar','mNu E0 Norm Bkg',... free parameter
    'NonPoissonScaleFactor',1.064,...
    'minuitOpt','min ; minos',...
    'FSDFlag','Sibille0p5eV',...
    'ELossFlag','KatrinT2',...
    'SysBudget',24,...
    'RadiativeFlag','ON',...
    'SynchrotronFlag','ON',...
    'AngularTFFlag','OFF',...
    'DopplerEffectFlag','OFF');
Real.exclDataStart = Real.GetexclDataStart(range);

%% stat + syst
Real.ResetFitResults;
Real.chi2 = 'chi2CMShape';
Real.NonPoissonScaleFactor  = 1.064;
Real.FitResult.par(1) = -0.8;
ScanResultsCM = Real.GetAsymFitError('Mode','Uniform','ParScanMax',1.8,'nFitMax',nFitMax,'SanityPlot','OFF');
%% stat. only
Real.ResetFitResults;
Real.chi2 = 'chi2Stat';
Real.FitResult.par(1) = -0.8;
Real.NonPoissonScaleFactor  = 1;
ScanResultsStat = Real.GetAsymFitError('Mode','Uniform','ParScanMax',1.8,'nFitMax',nFitMax,'SanityPlot','OFF');

%% plot
Real.chi2 = 'chi2CMShape';
pHandle = Real.PlotChi2Curve('Parameter','mNu','ScanResult',ScanResultsCM,'PlotBf','ON','HoldOn','OFF');
%Real.chi2 = 'chi2Stat';
%Real.PlotChi2Curve('Parameter','mNu','ScanResult',ScanResultsStat,'PlotBf','ON','HoldOn','ON');

%%
close all
mNuSq =[flipud(ScanResults.ParScan(:,2));ScanResults.ParScan(2:end,1)];
Chi2 = [flipud(ScanResults.chi2min(:,2));ScanResults.chi2min(2:end,1)];

plot(mNuSq,Chi2,':');