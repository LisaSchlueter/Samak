% Script to check impact of column density on endpoint:
RunList = 'FTpaper';
exclDataStart = 13;
WGTS_CD_MolPerCm2_V = 0.9:0.02:1.1;
RecomputeFlag = 'OFF';

switch exclDataStart
    case 1
        belowE0 = 1600;
    case 7 
        belowE0 = 400;
    case 9
        belowE0 = 200;
    case 12
        belowE0 = 125;
    case 13
        belowE0 = 100;
end
savedir = [getenv('SamakPath'),'First-Tritium-Paper/PlotScripts/results/'];
savename = sprintf('ColumnDensityImpactE0_%s_%.0feV_pm%.0fpercent_%.0fstep.mat',RunList,belowE0,...
    (max(WGTS_CD_MolPerCm2_V)-1)*100,(WGTS_CD_MolPerCm2_V(2)-WGTS_CD_MolPerCm2_V(1))*100);
saveresult = [savedir,savename];

if exist(saveresult,'file') && strcmp(RecomputeFlag,'OFF')
    load(saveresult);
else
     M = MultiRunAnalysis('RunList',RunList,'exclDataStart',exclDataStart,'ELossFlag','Abdurashitov',...
         'fixPar','1 5 6 7 8 9 10 11','FSDFlag','SAENZ','NonPoissonScaleFactor',1.0);
     [parRD, errRD, chi2RD, WGTS_CD_MolPerCm2_local, CD_bestfit] = ...
         M.RhoDScan('WGTS_CD_MolPerCm2_V',WGTS_CD_MolPerCm2_V,'saveplot','ON');
 save(saveresult,'parRD','errRD','chi2RD', 'WGTS_CD_MolPerCm2_local','CD_bestfit','WGTS_CD_MolPerCm2_V');
end

%% plot
f12 = figure(12);
set(f12, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.7]);

e2 = errorbar(WGTS_CD_MolPerCm2_local./WGTS_CD_MolPerCm2_local(WGTS_CD_MolPerCm2_V==1)*100,...
    parRD(2,:)+18573.7-mean(parRD(2,:)+18573.7),errRD(2,:),'o-', 'Color',rgb('IndianRed'),'LineWidth',2,...
    'MarkerFaceColor',rgb('IndianRed'),'MarkerSize',8);
xlabel('column density (%)');
ylabel('E_{0_{eff}} - < E_{0_{eff}}> (eV)');

rhomin = WGTS_CD_MolPerCm2_V(1)*100;
rhomax = WGTS_CD_MolPerCm2_V(end)*100;
xlim([rhomin rhomax]);
ylim([-0.65, 0.65])
xticks([rhomin:2:rhomax]);
ytickformat('%.5g')
grid on;
PrettyFigureFormat;
leg = legend(sprintf('100%% column density = %.3g mol/cm^2',...
    WGTS_CD_MolPerCm2_local(WGTS_CD_MolPerCm2_V==1)),'Location','northwest');
legend boxoff
set(gca,'FontSize',22);
print(gcf,[strrep(savedir,'results','plots'),strrep(savename,'mat','png')],'-dpng','-r450')
