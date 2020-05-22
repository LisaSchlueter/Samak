% Calculte and/or plot systematics breakdown for KNM1

% go into directory, where results are stored
mypath = [getenv('SamakPath'),'knm1ana/knm1sensitivity'];
system(['cd ',mypath]);
startup;
range = 40;
RunAnaArg = {'RunList','KNM1',...
    'fixPar','mNu E0 Norm Bkg',...
    'SysBudget',22,...
    'FSDFlag','Sibille0p5eV',...
    'fitter','minuit',...
    'minuitOpt','min;minos',...
    'NonPoissonScaleFactor',1.064,...
    'RadiativeFlag','ON',...
    'ELossFlag','KatrinT2',...
    'AngularTFFlag','OFF'};
%% Monte Carlo expectation
M = MultiRunAnalysis(RunAnaArg{:},'DataType','Twin');
M.exclDataStart=M.GetexclDataStart(range);
%%
S = RunSensitivity('RunAnaObj',M);
S.RecomputeFlag='OFF';
S.LimitFlag = 'Central';
S.ConfLevel=0; % 0 == 1 sigma

%S.PlotSysBreakdownBars('Ranges',14,'SavePlot','pdf','HoldOn','OFF','SysInfoBox','OFF','DispTitle','OFF');
S.PlotSysBreakdownBars2('Ranges',M.exclDataStart,'SavePlot','pdf','HoldOn','OFF','SysInfoBox','OFF');

 %% KNM1 real data
 D = MultiRunAnalysis(RunAnaArg{:},'DataType','Real');
Sd = RunSensitivity('RunAnaObj',D);
Sd.RecomputeFlag='OFF';
Sd.LimitFlag = 'Central';
Sd.ConfLevel=0; % 0 == 1 sigma
%Sd.PlotSysBreakdownBars('Ranges',14,'SavePlot','pdf','HoldOn','OFF','SysInfoBox','OFF','DispTitle','OFF');
Sd.PlotSysBreakdownBars2('Ranges',D.exclDataStart,'SavePlot','pdf','HoldOn','OFF','SysInfoBox','OFF');

% Sd.PlotSysBreakdownBars('Ranges',14,'SavePlot','ON','HoldOn','OFF','SysInfoBox','OFF','DispTitle','OFF');
% Sd.PlotSysBreakdownBars2('Ranges',14,'SavePlot','ON','HoldOn','OFF','SysInfoBox','OFF');
% 
% %structfun(@(x) round(x(1),2),Sd.MultiLpar,'UniformOutput',0)
% %structfun(@(x) (x(1).^2-S.MultiLpar.Stat(1)^2).^0.5,S.MultiLpar,'UniformOutput',0)
% %%  Comparision
% S.PlotCompareSensitivityTwinData;