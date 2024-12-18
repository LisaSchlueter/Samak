%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                       PLOT For Blessing                               %%
%%                                                                       %%
%%                      T. Lasserre, 18/03/2019                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run List
% 

%% 
SPR=MultiRunAnalysis(...
    'DataType','Real',...
    'RunList','KNM1_175mvRW',...
    'exclDataStart',1,...
    'fixPar','1 5 6 7 8 9 10',...
    'chi2','chi2Stat',...
    'ELossFlag','KatrinD2',...
    'StackTolerance',0.2,...
    'NonPoissonScaleFactor',1.5,...
    'Debug','ON',...
    'AnaFlag','StackPixel',...
    'PixList',[]);

%% Plot qU-distribution
%SPR.qUDistribution('saveplot','ON');

%% Plot Single-Run Informations
SPR.chi2='chi2Stat';
SPR.FitRunList('Recompute','ON')
SPR.FitRunList_SlowControl
SPR.PlotSCdistributions
SPR.PlotFitResultsDistributions

%% Systematics
SPR.chi2='chi2CMShape';
systematics = struct(...
    'RF_EL','ON',...
    'RF_BF','ON',...
    'RF_RX','ON',...
    'FSD','ON',...
    'TASR','ON',...
    'TCoff_RAD','ON',...
    'TCoff_OTHER','ON',...
    'DOPoff','OFF');
CMoptions   = {
    'WGTS_TASR_RelErr',0.001,...
    'FSDNorm_RelErr',0.01,...
    'FSDShapeGS_RelErr',0.04,...
    'FSDShapeES_RelErr',0.18,...
    'MACE_Ba_T_RelErr',0.01,...
    'MACE_Bmax_T_RelErr',0.002,...
    'WGTS_B_T_RelErr',0.025,...
    'WGTS_CD_MolPerCm2_RelErr',0.02,...
    'ISXsection_RelErr',0.02,...
    'BkgCM','ON',...
    'Stack','OFF'...
    };
SPR.ComputeCM('RecomputeFlag','OFF','SysEffects',systematics,CMoptions{:});

%% Fit
SPR.Fit('CATS','ON');
SPR.PlotFit('Mode','Rate','saveplot','png')

%% PlotDataModel
SPR.PlotDataModel_KNM1('saveplot','ON');
