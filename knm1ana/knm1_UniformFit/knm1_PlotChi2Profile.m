
% plot chi2 profiles
savedir = [getenv('SamakPath'),'tritium-data/fit/Knm1/Chi2Profile/Uniform/'];

% load CM
savenameCM = sprintf('%sChi2Profile_Real_UniformScan_mNu_Knm1_UniformFPD_chi2CMShape_SysBudget24_NP1.064_FitParE0BkgNorm_nFit50_SibilleFull_min-2.6_max1.mat',savedir);
dCM = importdata(savenameCM);
chi2_CM = zeros(99,1);
mNuSq_CM = zeros(99,1);
chi2_CM(1:50) = flipud(dCM.ScanResults.chi2min(:,2));
mNuSq_CM(1:50) = flipud(dCM.ScanResults.ParScan(:,2));
chi2_CM(51:99) = dCM.ScanResults.chi2min(2:end,1);
mNuSq_CM(51:99) = dCM.ScanResults.ParScan(2:end,1);

% load stat
savenameStat = sprintf('%sChi2Profile_Real_UniformScan_mNu_Knm1_UniformFPD_chi2Stat_NP1.000_FitParE0BkgNorm_nFit50_SibilleFull_min-2.6_max1.mat',savedir);
dStat = importdata(savenameStat);
chi2_Stat = zeros(99,1);
mNuSq_Stat = zeros(99,1);
chi2_Stat(1:50) = flipud(dStat.ScanResults.chi2min(:,2));
mNuSq_Stat(1:50) = flipud(dStat.ScanResults.ParScan(:,2));
chi2_Stat(51:99) = dStat.ScanResults.chi2min(2:end,1);
mNuSq_Stat(51:99) = dStat.ScanResults.ParScan(2:end,1);

%% plot

close all
figureHandle = figure('Units','normalized','Position',[0.1,0.1,0.5,0.45]);
%GetFigure
acm = area([dCM.ParScanDown dCM.ParScanUp],[50,50],'FaceColor',rgb('DodgerBlue'),'FaceAlpha',0.3,'EdgeColor',rgb('DodgerBlue'));
hold on;
astat = area([dStat.ParScanDown dStat.ParScanUp],[50,50],'FaceColor',rgb('Orange'),'FaceAlpha',0.2,'EdgeColor',rgb('Orange'));

pcm = plot(mNuSq_CM,chi2_CM,'-','Color',rgb('DodgerBlue'),'LineWidth',3);
pstat = plot(mNuSq_Stat,chi2_Stat,'-.','Color',rgb('Orange'),'LineWidth',3);
plot(dStat.ScanResults.BestFit.par(1),dStat.ScanResults.BestFit.chi2,'.','MarkerSize',20,'Color',rgb('Orange'));
plot(dCM.ScanResults.BestFit.par(1),dCM.ScanResults.BestFit.chi2,'.','MarkerSize',20,'Color',rgb('DodgerBlue'));
ylim([21.2 25.2]);
xlim([-2.6 0.6]);
PrettyFigureFormat('FontSize',22);
xlabel(sprintf('{\\itm}_\\nu^2 (eV^2)'));
ylabel(sprintf('\\chi^2 (%.0f dof)',23));

leg = legend([pstat,pcm],'Stat.','Stat. and syst.','Location','north');
PrettyLegendFormat(leg);
leg.NumColumns =2;
leg.FontSize = get(gca,'FontSize');

%% save plot dir
pltdir =  [getenv('SamakPath'),'knm1ana/knm1_UniformFit/plots/'];
MakeDir(pltdir);
pltname =sprintf('%sknm1_PlotChi2Profile.png',pltdir);
print(pltname,'-dpng','-r400');



