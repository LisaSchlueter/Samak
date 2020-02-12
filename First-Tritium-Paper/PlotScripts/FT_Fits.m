% script to fit and plot fits for FT paper
FT = MultiRunAnalysis('RunList','FTpaper','chi2','chi2CMShape',...
    'ELossFlag','Abdurashitov','exclDataStart',7,'SysBudget',0,...%0
    'DataType','Real','FSDFlag','SAENZ',...
    'fixPar','1 5 6 7 8 9 10 11','NonPoissonScaleFactor',1,...
    'minuitOpt','min;migrad');
Q_i = FT.ModelObj.Q_i;
%% Fit
exclDataStartAll = [13,9,7]; % fit ranges: 1600, 400,200,125
for i=1:numel(exclDataStartAll)
    FT.exclDataStart=exclDataStartAll(i);
    FT.ComputeCM;
    FT.Fit
    FT.ModelObj.Q_i = 18574;
    FT.PlotFit('SavePlot','pdf','YLimRes',[-2.5 2.5],'ErrorBarScaling',50,'FitResultsFlag','OFF');
    FT.ModelObj.Q_i = Q_i;
end
