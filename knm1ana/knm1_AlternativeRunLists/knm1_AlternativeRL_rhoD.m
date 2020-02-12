%
% Determine Alternative Run Lists for checking
% the KNM1 Stability of the Analysis wrt RhoD
%
% Thierry Lasserre
% Last updated: 11/07/2019
%

%%
[Real , Twin] = Knm1RealTwin_Create();

%% Plot RhoD TimeLine
Real.PlotSCdata_RhoD;

%% WGTS_CD_MolPerCm2 Distribution
RhoD_StackedRuns    = Real.SingleRunData.WGTS_CD_MolPerCm2;

%% Alternative Run List
m = mean(RhoD_StackedRuns);
s = std(RhoD_StackedRuns);
Real.SingleRunData.Select_all = (RhoD_StackedRuns<m+1*s & RhoD_StackedRuns>m-1*s);
RhoD_StackedRuns    = Real.SingleRunData.WGTS_CD_MolPerCm2(Real.SingleRunData.Select_all);
RhoD_StackedRunsErr = std(Real.SingleRunData.WGTS_CD_MolPerCm2_SubRun(:,Real.SingleRunData.Select_all));
RhoD_NotStackedRuns = Real.SingleRunData.WGTS_CD_MolPerCm2(:,Real.SingleRunData.Select_all==0);

%% Plot Data
fig10100 = figure(10100);
set(fig10100, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.9, 0.6]);
subplot(1,4,[1 3])
hstack=plot(Real.SingleRunData.StartTimeStamp(Real.SingleRunData.Select_all)',...
    RhoD_StackedRuns,'ks-.','MarkerSize',8,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',1);
hold on
hold off

xlabel('Time');
ylabel('\rho d (mol/cm^2)');
legend([hstack],sprintf('<\\rho d> = %.3g mol/cm^2',Real.RunData.WGTS_CD_MolPerCm2));
legend boxoff;
title('Column Density from WGTS Throughput');
%datetick('x', 'dd-mm-yyyy');
xtickangle(45);
PrettyFigureFormat;set(gca,'FontSize',16);

subplot(1,4,4)
fiths = histfit(RhoD_StackedRuns,5);
fiths(1).FaceColor = rgb('IndianRed');
fiths(1).LineStyle = '-';
xlabel('\rho d');
xtickangle(45);
ylabel('Runs');
legend([fiths],sprintf('%.0f runs',numel(Real.RunList(Real.SingleRunData.Select_all))),'Location','NorthWest'); legend boxoff;
PrettyFigureFormat;set(gca,'FontSize',16);

%% Run Lists
fprintf('KNM1_RhoD Run List:\n');
KNM1_RhoD     = Real.RunList(Real.SingleRunData.Select_all);
disp(num2str(KNM1_RhoD'));
fprintf('KNM1_AntiRhoD Run List:\n');
KNM1_antiRhoD = Real.RunList(~Real.SingleRunData.Select_all);
disp(num2str(KNM1_antiRhoD'));

% Save
savefile=sprintf('plots/MRA_SC%0.fdata_RhoD_AleternativeRL.png',numel(Real.RunList(Real.SingleRunData.Select_all)));
export_fig(gcf,savefile,'-q101','-m3');

%% Load Analysis
options = {...
    'RunList','KNM1_RhoD',...
    'exclDataStart',14,...
    'fixPar','5 6 7 8 9 10 11',...
    'chi2','chi2Stat',...
    'ELossFlag','KatrinD2',...
    'StackTolerance',1,...
    'NonPoissonScaleFactor',1.1,...
    'Debug','ON',...
    'AnaFlag','StackPixel',...
    'PixList',[]};
KNM1_RhoD = MultiRunAnalysis('DataType','Twin',options{:});

if KNM1_RhoD.DataType=='Real'
    KNM1_RhoD.PlotSCdistributions;
    KNM1_RhoD.PlotSCdata_RhoD;
end

%% Sensitivity
if KNM1_RhoD.DataType=='Twin'
KNM1_RhoD.Fit;KNM1_RhoD.PlotFit;
KNM1_RhoD.PlotDataModel_KNM1;
end