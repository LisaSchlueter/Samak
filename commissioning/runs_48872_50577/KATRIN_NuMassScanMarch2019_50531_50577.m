%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                       PLOT For Blessing                               %%
%%                                                                       %%
%%                      T. Lasserre, 18/03/2019                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run List
% 50531 ... 50577

%% 
SPR=MultiRunAnalysis(...
    'DataType','Real',...
    'RunList','March2019',...
    'exclDataStart',2,...
    'fixPar','1 5 6 7 8 9 10',...
    'chi2','chi2Stat',...
    'ELossFlag','KatrinD2',...
    'StackTolerance',0.2,...
    'NonPoissonScaleFactor',1.5,...
    'Debug','ON');

%% Plot qU-distribution
SPR.qUDistribution('saveplot','ON');

% Plot Single-Run Informations
SPR.chi2='chi2Stat';
SPR.FitRunList('Recompute','ON')
SPR.FitRunList_SlowControl
SPR.PlotFitResultsDistributions
