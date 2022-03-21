
% plot multi-ring chi2 profiles
savedir = [getenv('SamakPath'),'tritium-data/fit/Knm2/Chi2Profile/Ring_Full/'];

% load CM
savenameCM = sprintf('Chi2Profile_Real_UniformScan_mNu_Knm2_Ring_FullFPD_chi2CMShape_SysBudget41_NP1.112_FitParE0BkgNormqU_nFit20_min-2_max1.mat',savedir);
dCM = importdata(savenameCM);
nFits = numel(dCM.ScanResults.chi2min);
chi2_CM = zeros(nFits-1,1);
mNuSq_CM = zeros(nFits-1,1);
chi2_CM(1:nFits/2) = flipud(dCM.ScanResults.chi2min(:,2));
mNuSq_CM(1:nFits/2) = flipud(dCM.ScanResults.ParScan(:,2));
chi2_CM(nFits/2+1:nFits-1) = dCM.ScanResults.chi2min(2:end,1);
mNuSq_CM(nFits/2+1:nFits-1) = dCM.ScanResults.ParScan(2:end,1);

%load stat
savenameStat = sprintf('%sChi2Profile_Real_Parallel_mNu_Knm2_Ring_FullFPD_chi2Stat_NP1.000_FitParmNuE0BkgNormqU_nFit30_Min-0.5_Max1_KNM2_0p1eV_40eVrange.mat',savedir);
dStat = importdata(savenameStat);
chi2_Stat = dStat.chi2min;
mNuSq_Stat = dStat.ParScan;


%% get 1 sigma uncertainties CM
 mNuSq_CM_inter = linspace(min(mNuSq_CM),max(mNuSq_CM),5e3);
 chi2_CM_inter = interp1(mNuSq_CM,chi2_CM,mNuSq_CM_inter,'spline');
 chi2_CM_bf = min(chi2_CM_inter);
  Idx_bf = find(chi2_CM_inter== chi2_CM_bf);
 mNuSq_CM_bf = mNuSq_CM_inter(Idx_bf);

 mNuSq_CM_low = interp1(chi2_CM_inter(1:Idx_bf),mNuSq_CM_inter(1:Idx_bf), chi2_CM_bf+1,'spline');
 mNuSq_CM_up = interp1(chi2_CM_inter(Idx_bf:end),mNuSq_CM_inter(Idx_bf:end), chi2_CM_bf+1,'spline');
 
%% plot

close all
figureHandle = figure('Units','normalized','Position',[0.1,0.1,0.5,0.45]);
%GetFigure
acm = area([mNuSq_CM_low mNuSq_CM_up],[150,150],'FaceColor',rgb('DodgerBlue'),'FaceAlpha',0.3,'EdgeColor',rgb('DodgerBlue'));
hold on;
astat = area([dStat.BestFit.parLowLim dStat.BestFit.parUpLim],[150,150],'FaceColor',rgb('Orange'),'FaceAlpha',0.2,'EdgeColor',rgb('Orange'));

pcm = plot(mNuSq_CM_inter,chi2_CM_inter,'-','Color',rgb('DodgerBlue'),'LineWidth',3);
pstat = plot(mNuSq_Stat,chi2_Stat,'-.','Color',rgb('Orange'),'LineWidth',3);

plot(dStat.BestFit.par,dStat.BestFit.chi2,'.','MarkerSize',20,'Color',rgb('Orange'));
plot(mNuSq_CM_bf, chi2_CM_bf,'.','MarkerSize',20,'Color',rgb('DodgerBlue'));
ylim([min(chi2_Stat)-1 min(chi2_CM)+5.4]);
xlim(mNuSq_CM_bf+[-1 1].*0.737);
PrettyFigureFormat('FontSize',22);
xlabel(sprintf('{\\itm}_\\nu^{ 2} (eV^{ 2})'));
ylabel(sprintf('\\chi^2 (%.0f dof)',dStat.dof(1)-1));

leg = legend([pstat,pcm],'Stat.','Stat. and syst.','Location','north');
PrettyLegendFormat(leg);
leg.NumColumns =2;
leg.FontSize = get(gca,'FontSize');

%% save plot dir
pltdir =  [getenv('SamakPath'),'knm2ana/knm2_PngBkg/plots/'];
MakeDir(pltdir);
pltname =sprintf('%sknm2_PlotChi2Profile.pdf',pltdir);
export_fig(pltname);
%print(pltname,'-dpng','-r400');



