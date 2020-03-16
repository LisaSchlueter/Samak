% Calculte and/or plot systematics breakdown for KNM2

range = 40;
% go into directory, where results are stored
mypath = [getenv('SamakPath'),'knm2ana/knm2_sensitivity'];
system(['cd ',mypath]);
startup;

RunAnaArg = {'RunList','KNM2_Prompt',...
    'fixPar','mNu E0 Norm Bkg',...
    'SysBudget',33,... % 31= knm2 preliminary input
    'FSDFlag','BlindingKNM2',...
    'fitter','matlab',...
    'NonPoissonScaleFactor',1.0877,...
    'RadiativeFlag','ON',...
    'ELossFlag','KatrinT2'};
%% Monte Carlo expectation
TwinOpt = {'DataType','Twin','TwinBias_Q',18573.7,'TwinBias_mnuSq',0};
MC = MultiRunAnalysis(RunAnaArg{:},TwinOpt{:});
MC.exclDataStart = MC.GetexclDataStart(range);
%%
S = RunSensitivity('RunAnaObj',MC);
S.RecomputeFlag='ON';
S.LimitFlag = 'Central';
S.ConfLevel=0; % 0 == 1 sigma

%S.PlotSysBreakdownBars('Ranges',14,'SavePlot','pdf','HoldOn','OFF','SysInfoBox','OFF','DispTitle','OFF');
S.PlotSysBreakdownBars2('Ranges',11,'SavePlot','pdf','HoldOn','OFF','SysInfoBox','OFF');


 %% KNM1 real data
% D = MultiRunAnalysis(RunAnaArg{:},'DataType','Real');
% Sd = RunSensitivity('RunAnaObj',D);
% Sd.RecomputeFlag='OFF';
% Sd.LimitFlag = 'Central';
% Sd.ConfLevel=0; % 0 == 1 sigma
% Sd.PlotSysBreakdownBars('Ranges',14,'SavePlot','pdf','HoldOn','OFF','SysInfoBox','OFF','DispTitle','OFF');
% Sd.PlotSysBreakdownBars2('Ranges',14,'SavePlot','pdf','HoldOn','OFF','SysInfoBox','OFF');
% 
% Sd.PlotSysBreakdownBars('Ranges',14,'SavePlot','ON','HoldOn','OFF','SysInfoBox','OFF','DispTitle','OFF');
% Sd.PlotSysBreakdownBars2('Ranges',14,'SavePlot','ON','HoldOn','OFF','SysInfoBox','OFF');
% 
% %structfun(@(x) round(x(1),2),Sd.MultiLpar,'UniformOutput',0)
% %structfun(@(x) (x(1).^2-S.MultiLpar.Stat(1)^2).^0.5,S.MultiLpar,'UniformOutput',0)
% %%  Comparision
% S.PlotCompareSensitivityTwinData;