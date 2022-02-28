% twins for the three periods seperatetly 

range = 40;
BKG_PtSlope = 3*1e-06;

RunAnaArg = {'RunList','KNM2_Prompt',...
    'fixPar','E0 Bkg Norm',...
    'DataType','Real',...  
    'FSDFlag','KNM2_0p1eV',...
    'ELossFlag','KatrinT2A20',...
    'AnaFlag','StackPixel',...
    'BKG_PtSlope',3*1e-06,...
    'chi2','chi2Stat',...
    'DopplerEffectFlag','FSD',...
    'RadiativeFlag','ON',...
    'NonPoissonScaleFactor',1,...
    'FSD_Sigma',sqrt(0.0124+0.0025),...
    'BKG_PtSlope',BKG_PtSlope};

% read data and set up model
M = MultiRunAnalysis(RunAnaArg{:});
M.exclDataStart = M.GetexclDataStart(range);
%% fit all runs one by one (or load fit results from file)
M.FitRunList;

return
%%
M.PlotFitRunList('DisplayStyle','Abs','Parameter','RhoD','HideGaps','OFF','YLim',[4.20,4.25].*1e17)
M.PlotFitRunList('DisplayStyle','Abs','Parameter','T2','HideGaps','OFF','ShowRWPeriods','ON','Saveplot','pdf');
%% plots
M.PlotFitRunList('saveplot','pdf','DisplayStyle','Rel','Parameter','E0','YLim',[-0.65,0.8],'DispHist','ON');
M.PlotFitRunList('saveplot','pdf','DisplayStyle','Abs','Parameter','B','YLim','','DispHist','ON');
M.PlotFitRunList('saveplot','pdf','DisplayStyle','Abs','Parameter','N','YLim',[0.89,1.05],'DispHist','ON');