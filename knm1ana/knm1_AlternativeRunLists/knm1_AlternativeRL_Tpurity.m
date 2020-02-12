%
% Determine Alternative Run Lists for checking
% the KNM1 Stability of the Analysis wrt Tritium Purity
%
% Thierry Lasserre
% Last updated: 1Z/07/2019
%

%%
[Real , Twin] = Knm1RealTwin_Create();

%% Plot RhoD TimeLine
Real.PlotSCdata_TritiumPurity;

%% Select Runs
epsT_StackedRuns    = Real.SingleRunData.WGTS_MolFrac_TT+0.5*Real.SingleRunData.WGTS_MolFrac_HT+0.5*Real.SingleRunData.WGTS_MolFrac_DT;
m                             = median(epsT_StackedRuns);
Real.SingleRunData.Select_all = epsT_StackedRuns>m;
epsT_StackedRuns              = epsT_StackedRuns(Real.SingleRunData.Select_all);

%% Plot Data
fig10100 = figure(10100);
set(fig10100, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.9, 0.6]);
subplot(1,4,[1 3])
hstack=plot(Real.SingleRunData.StartTimeStamp(Real.SingleRunData.Select_all)',...
    epsT_StackedRuns,'ks-.','MarkerSize',8,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',1);
hold on
hold off

xlabel('Time');
ylabel('\rho d (mol/cm^2)');
legend([hstack],sprintf('<tritium purity> = %.3g',Real.RunData.WGTS_MolFrac_TT+0.5*Real.RunData.WGTS_MolFrac_HT+0.5*Real.RunData.WGTS_MolFrac_DT));
legend boxoff;
title('Tritium Purity');
%datetick('x', 'dd-mm-yyyy');
xtickangle(45);
PrettyFigureFormat;set(gca,'FontSize',16);

subplot(1,4,4)
fiths = histfit(epsT_StackedRuns,5);
fiths(1).FaceColor = rgb('IndianRed');
fiths(1).LineStyle = '-';
xlabel('\rho d');
xtickangle(45);
ylabel('Runs');
legend([fiths],sprintf('%.0f runs',numel(Real.RunList(Real.SingleRunData.Select_all))),'Location','NorthWest'); legend boxoff;
PrettyFigureFormat;set(gca,'FontSize',16);

% Run Lists
fprintf('KNM1-TpurityHigh Run List:\n');
KNM1TpurityHigh     = Real.RunList(Real.SingleRunData.Select_all);
disp(num2str(KNM1TpurityHigh'));
fprintf('KNM1-TpurityLow Run List:\n');
KNM1TpurityLow      = Real.RunList(~Real.SingleRunData.Select_all);
disp(num2str(KNM1TpurityLow'));

% Save
savefile=sprintf('plots/MRA_SC%0.fdata_Tpurity_AleternativeRL.png',numel(Real.RunList(Real.SingleRunData.Select_all)));
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