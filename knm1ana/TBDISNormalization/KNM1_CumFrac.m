% Create MultiRunAnalysis Object
if ~exist('SPR_R','var') || ~exist('SPR_T','var')
    [SPR_R , SPR_T] = Knm1RealTwin_Create();
end

% Get Single Run Model Objects
if isempty(SPR_R.SingleRunObj)
SPR_R.LoadSingleRunObj;
end

% Extract Cum Frac
SR_CumFrac = cellfun(@(x) x.CumFrac,SPR_R.SingleRunObj,'UniformOutput', true);

% Plot
CumFracFig = figure(1);
set(CumFracFig,'Units', 'normalized','Position', [0.1, 0.1, 0.8, 0.8]);
s1=subplot(2,1,1)
title('Fraction of Decays in TBD Tail (CumFrac)')
stairs(SR_CumFrac,'LineWidth',2);
xlabel('Run');
ylabel('fraction');
PrettyFigureFormat
title('Fraction of Decays in TBD Tail (CumFrac)')
s2=subplot(2,1,2)
stairs(SR_CumFrac./mean(SR_CumFrac),'LineWidth',2);
xlabel('Run');
ylabel('fraction/mean(fraction)');
PrettyFigureFormat
export_fig(1,'plots/knm1_cumfrac_1.png','-r300');

%% Plot  CumFrac VS qUmin
CumFracFig = figure(2);
set(CumFracFig,'Units', 'normalized','Position', [0.1, 0.1, 0.8, 0.8]);
scatter(min(SPR_R.SingleRunData.qU),SR_CumFrac,'s','LineWidth',2);
ylabel('CumFrac');
xlabel('min(qU)');
title('Fraction of Decays in TBD Tail (CumFrac) Versus min(qU)');
PrettyFigureFormat
export_fig(2,'plots/knm1_cumfrac_2.png','-r300');

%% Plot  CumFrac VS RhoD
CumFracFig = figure(3);
set(CumFracFig,'Units', 'normalized','Position', [0.1, 0.1, 0.8, 0.8]);
scatter(SPR_R.SingleRunData.WGTS_CD_MolPerCm2,SR_CumFrac','s','LineWidth',2);
ylabel('CumFrac');
xlabel('\rho d');
title('Fraction of Decays in TBD Tail (CumFrac) Versus Column Density');
PrettyFigureFormat
export_fig(3,'plots/knm1_cumfrac_3.png','-r300');

%% Plot  CumFrac VS [TT] / [HT] / [DT]
CumFracFig = figure(4);
set(CumFracFig,'Units', 'normalized','Position', [0.1, 0.1, 0.8, 0.8]);
s1=subplot(3,1,1)
scatter(SPR_R.SingleRunData.WGTS_MolFrac_TT,SR_CumFrac','s','LineWidth',2);
ylabel('CumFrac');
xlabel('[TT]');
PrettyFigureFormat
s2=subplot(3,1,2)
scatter(SPR_R.SingleRunData.WGTS_MolFrac_HT,SR_CumFrac','s','LineWidth',2);
ylabel('CumFrac');
xlabel('[HT]');
PrettyFigureFormat
s2=subplot(3,1,3)
scatter(SPR_R.SingleRunData.WGTS_MolFrac_DT,SR_CumFrac','s','LineWidth',2);
ylabel('CumFrac');
xlabel('[DT]');
PrettyFigureFormat
export_fig(4,'plots/knm1_cumfrac_4.png','-r300');

%% Values
fprintf('Mean CumFrac = %g\n',mean(SR_CumFrac));
[~ , i ] = min(SR_CumFrac);
fprintf('Minimum CumFrac = %g for Run %.0f \n',min(SR_CumFrac), SPR_R.RunList(i));
[~ , i ] = max(SR_CumFrac);
fprintf('Minimum CumFrac = %g for Run %.0f \n',max(SR_CumFrac), SPR_R.RunList(i));

%% CumFrac Versus Fit Parameters
CumFracFig = figure(5);
set(CumFracFig,'Units', 'normalized','Position', [0.1, 0.1, 0.8, 0.8]);
subplot(2,2,1)
scatter(SPR_R.SingleRun_FitResults.chi2Stat.N,SR_CumFrac','s','LineWidth',2);
ylabel('CumFrac');
xlabel('Single Run Fit Normalization');
PrettyFigureFormat
subplot(2,2,2)
scatter(SPR_R.SingleRun_FitResults.chi2Stat.E0,SR_CumFrac','s','LineWidth',2);
ylabel('CumFrac');
xlabel('Single Run Fit E0');
PrettyFigureFormat
subplot(2,2,3)
scatter(SPR_R.SingleRun_FitResults.chi2Stat.B,SR_CumFrac','s','LineWidth',2);
ylabel('CumFrac');
xlabel('Single Run Fit Background (cps)');
PrettyFigureFormat
subplot(2,2,4)
scatter(SPR_R.SingleRun_FitResults.chi2Stat.chi2min,SR_CumFrac','s','LineWidth',2);
ylabel('CumFrac');
xlabel('Single Run Fit \chi ^2');
PrettyFigureFormat
export_fig(5,'plots/knm1_cumfrac_5.png','-r300');

