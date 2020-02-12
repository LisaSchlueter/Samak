% twins for the three periods seperatetly 
RunList = 'KNM2_RW1'; % 
fixPar = 'E0 Bkg Norm'; % free parameter
DataType = 'Real';
FSDFlag = 'BlindingKNM2';
ELossFlag = 'KatrinT2';
AnaFlag = 'StackPixel'; % uniform FPD
exclDataStart = 11;
chi2 = 'chi2Stat';

RunAnaArg = {'RunList',RunList,...
    'fixPar',fixPar,...
    'DataType',DataType,...
    'FSDFlag',FSDFlag,...
    'ELossFlag',ELossFlag,...
    'AnaFlag',AnaFlag,'chi2',chi2,...
    'TwinBias_Q','Fit',...
    'exclDataStart',exclDataStart};

% read data and set up model
M = MultiRunAnalysis(RunAnaArg{:});

%% fit all runs one by one (or load fit results from file)
M.FitRunList;

%% plots
M.PlotFitRunList('saveplot','pdf','DisplayStyle','Rel','Parameter','E0','YLim',[-0.65,0.8],'DispHist','ON');
M.PlotFitRunList('saveplot','pdf','DisplayStyle','Abs','Parameter','B','YLim','','DispHist','ON');
M.PlotFitRunList('saveplot','pdf','DisplayStyle','Abs','Parameter','N','YLim',[0.89,1.05],'DispHist','ON');