RunList = 'KNM1';
DataType = 'Real';
chi2 = 'chi2Stat';
fixPar = '1 5 6 7 8 9 10 11';
exclDataStart = 14;

FSDFlag ='SibilleFull';
M = MultiRunAnalysis('RunList',RunList,...
    'chi2','chi2Stat','DataType',DataType,...
    'fixPar',fixPar,...
    'RadiativeFlag','ON',...
    'NonPoissonScaleFactor',1.064,...
    'minuitOpt','min ; migrad',...
    'FSDFlag',FSDFlag,...
    'ELossFlag','KatrinT2',...
    'SysBudget',22,...
    'exclDataStart',exclDataStart);
%%
M.PlotFitRunList('SavePlot','ON','Parameter','E0','DisplayStyle','Rel','DispHist','ON','YLim',[-1,1.2]);
M.PlotFitRunList('SavePlot','ON','Parameter','B','DisplayStyle','Rel','DispHist','ON','YLim',[-40,55]);