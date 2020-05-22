% Test SterileAnalysis class
% Lisa, May2020
% plot qU Scan
% fixed parameter
%% settings for runanalysis
DataType = 'Twin';
%%
RunAnaArg = {'RunList','KNM1',...
    'fixPar','E0 Norm Bkg',...
    'DataType',DataType,...
    'FSDFlag','SibilleFull',...
    'ELossFlag','KatrinT2',...
    'AnaFlag','StackPixel',...
    'chi2','chi2Stat',...
    'ROIFlag','Default',...
    'SynchrotronFlag','ON',...
    'AngularTFFlag','OFF',...
    'ISCSFlag','Edep',...
    'TwinBias_Q',18573.73,...
    'SysBudget',24,...
    'pullFlag',99,...
    'NonPoissonScaleFactor',1};

T = MultiRunAnalysis(RunAnaArg{:});
T.chi2 = 'chi2CMShape';
%% settings sterile class
SterileArg = {'RunAnaObj',T,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
    'nGridSteps',50,...
    'SmartGrid','OFF',...
    'RecomputeFlag','OFF',...
    'SysEffect','all',...
    'RandMC','OFF',...
    'range',95,...
    'ConfLevel',95};

S = SterileAnalysis(SterileArg{:});
%% Samak plots
S.InterpMode = 'spline';
S.RunAnaObj.chi2 = 'chi2CMShape';
S.SysEffect = 'FSD';
S.nGridSteps = 25;

S.RunAnaObj.SysBudget = 241; % only FSD onset
S.LoadGridFile('CheckSmallerN','ON');
S.Interp1Grid('RecomputeFlag','ON');
pFSDonset = S.ContourPlot('CL',S.ConfLevel,'HoldOn','OFF',...
    'Color',rgb('DodgerBlue'),'LineStyle','-','BestFit','OFF','PlotSplines','OFF');

S.RunAnaObj.SysBudget = 242; % only FSD bin-to-bin
S.LoadGridFile('CheckSmallerN','ON');
S.Interp1Grid('RecomputeFlag','ON');
pFSDuncorr = S.ContourPlot('CL',S.ConfLevel,'HoldOn','ON',...
    'Color',rgb('FireBrick'),'LineStyle','-','BestFit','OFF','PlotSplines','OFF');

% S.SysEffect = 'FSD';
% S.RunAnaObj.SysBudget = 24; %
% S.LoadGridFile('CheckSmallerN','ON');
% S.Interp1Grid('RecomputeFlag','ON');
% pSysall = S.ContourPlot('CL',S.ConfLevel,'HoldOn','ON',...
%     'Color',rgb('DarkSlateGray'),'LineStyle','-','BestFit','OFF','PlotSplines','OFF');

S.SysEffect = 'all';
S.RunAnaObj.SysBudget = 24; %
S.LoadGridFile('CheckSmallerN','ON');
S.Interp1Grid('RecomputeFlag','ON');
pSysall = S.ContourPlot('CL',S.ConfLevel,'HoldOn','ON',...
    'Color',rgb('DimGray'),'LineStyle','-','BestFit','OFF','PlotSplines','OFF');

% S.RunAnaObj.chi2 = 'chi2Stat';
% S.LoadGridFile('CheckSmallerN','ON');
% S.Interp1Grid('RecomputeFlag','ON');
% pStat = S.ContourPlot('CL',S.ConfLevel,'HoldOn','ON',...
%     'Color',rgb('DimGray'),'LineStyle',':','BestFit','OFF','PlotSplines','OFF');

% fitrium
LineWidth = 2.5;
savedirF = [getenv('SamakPath'),'SterileAnalysis/GridSearchFiles/Knm1/Others/'];
fonset = sprintf('%scontour_KSN1_Fitrium_%s_%.0feV_FSD_onsetOnly_95CL_0.txt',savedirF,S.RunAnaObj.DataType,S.range);
dfFSDonset = importdata(fonset);

funcorr     = sprintf('%scontour_KSN1_Fitrium_%s_%.0feV_FSD_binsOnly_95CL_0.txt',savedirF,S.RunAnaObj.DataType,S.range);
dfFSDuncorr = importdata(funcorr);

ftotal     = sprintf('%scontour_KSN1_Fitrium_%s_%.0feV_total_95CL_0.txt',savedirF,S.RunAnaObj.DataType,S.range);
dftotal = importdata(ftotal);
pFtotal = plot(dftotal.data(:,1),dftotal.data(:,2),'LineStyle','-.','Color',rgb('Silver'),'LineWidth',LineWidth);
% % 
% fstat     = sprintf('%scontour_KSN1_Fitrium_%s_%.0feV_stat_95CL_0.txt',savedirF,S.RunAnaObj.DataType,S.range);
% dfstat = importdata(fstat);
% pFstat = plot(dfstat.data(:,1),dfstat.data(:,2),'LineStyle','-.','Color',rgb('Silver'),'LineWidth',LineWidth);
% % 
pFFSDonset = plot(dfFSDonset.data(:,1),dfFSDonset.data(:,2),'LineStyle','-.','Color',rgb('PowderBlue'),'LineWidth',LineWidth);
hold on;
pFFSDuncorr = plot(dfFSDuncorr.data(:,1),dfFSDuncorr.data(:,2),'LineStyle','-.','Color',rgb('Orange'),'LineWidth',LineWidth);

% %% 
legend([pFSDonset,pFFSDonset,pFSDuncorr,pFFSDuncorr],...
    'Samak (stat + FSD onset)', 'Fitrium (stat + FSD onset)',...
    'Samak (stat + FSD shape)', 'Fitrium (stat + FSD shape)',...
    'EdgeColor',rgb('Silver'),'Location','southwest');

%
% legend([pFSDonset,pFSDuncorr,pSysall,pFtotal],...
%     'Samak (stat + FSD onset)',...
%     'Samak (stat + FSD shape)', 'Samak (total)',...
%     'Fitrium (total)',...
%     'EdgeColor',rgb('Silver'),'Location','southwest');
xlim([1e-03,0.4])
ylim([3 90^2]);

title(sprintf('%s - %.0f eV range - %.0f%% C.L.',S.GetPlotTitle('Mode','data'),S.range,S.ConfLevel),...
    'FontWeight','normal','FontSize',get(gca,'FontSize'));

print(gcf,'KSN1_FSDuncertainty95eV_95CL.png','-dpng','-r400')