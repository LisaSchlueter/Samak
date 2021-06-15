RasterScan = 'ON';
DataType = 'Twin';
nGridSteps = 30;
range = 40;
savedir = [getenv('SamakPath'),'ksn2ana/ksn2_Systematics/results/'];


savename = sprintf('%sksn2_SystBreakdown_StatOverSyst_%s_%.0feV_RasterScan',savedir,DataType,range);

savenameRaster = [savename,'ON.mat'];
dRaster  = importdata(savenameRaster);
fprintf('load file %s \n',savenameRaster);
savenameScan = [savename,'OFF.mat'];
dScan  = importdata(savenameScan);
fprintf('load file %s \n',savenameScan);

%% plot both
GetFigure;
pRaster = plot(dRaster.mNu4Sq(1,:),dRaster.Ratio(1,:),'LineWidth',3,'Color',rgb('DodgerBlue'),'LineStyle','-');
hold on;
pScan = plot(dScan.mNu4Sq(1,:),dScan.Ratio(1,:),'LineWidth',3,'Color',rgb('IndianRed'),'LineStyle','-.');
ylabel(sprintf('\\sigma_{syst.}^2 / \\sigma_{total}^2'))
xlabel(sprintf('{\\itm}_4^2 (eV^2)'));
set(gca,'XScale','log');
PrettyFigureFormat('FontSize',22);
xlim([1 40^2]);

leg = legend([pRaster,pScan],'Raster scan (1 dof)','GridSearch (2 dof)');
PrettyLegendFormat(leg);
leg.Location = 'northwest';

plotdir = strrep(savedir,'results','plots');
plotname = sprintf('%sksn2_SystRatio_CompareRasterScanGridSearch.png',plotdir);
print(gcf,plotname,'-dpng','-r400');
fprintf('save plot to %s \n',plotname);
%% contour
GetFigure;
pRaster     = plot(dRaster.sin2t4_Tot(1,:),dRaster.mNu4Sq(1,:),'LineWidth',3,'Color',rgb('DodgerBlue'),'LineStyle','-');
hold on;
%pRasterStat = plot(dRaster.sin2t4_Stat(1,:),dRaster.mNu4Sq(1,:),'LineWidth',2,'Color',rgb('SkyBlue'),'LineStyle','-');

pScan = plot(dScan.sin2t4_Tot(1,:),dScan.mNu4Sq(1,:),'LineWidth',3,'Color',rgb('IndianRed'),'LineStyle','-.');
%pScanStat = plot(dScan.sin2t4_Stat(1,:),dScan.mNu4Sq(1,:),'LineWidth',2,'Color',rgb('Orange'),'LineStyle','-.');

xlabel(sprintf('|{\\itU}_{e4}|^2 '))
ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
set(gca,'XScale','log');
PrettyFigureFormat('FontSize',22);
ylim([1 40^2]);
xlim([3e-03 0.5]);
leg = legend([pRaster,pScan],sprintf('Raster scan (\\chi^2_{crit.} =  %.2f at 95%% C.L.)',GetDeltaChi2(95,1)),...
                             sprintf('Grid search (\\chi^2_{crit.} =  %.2f at 95%% C.L.)',GetDeltaChi2(95,2)));
% leg = legend([pRaster,pScan,pRasterStat,pScanStat],'Raster scan (total)','Grid search (total)',...
%                                                 'Raster scan (stat. only)','Grid search (stat. only)');
                                            
PrettyLegendFormat(leg);
leg.Location = 'southwest';
set(gca,'XScale','log');
set(gca,'YScale','log');
title(sprintf('Twin 40 eV range {\\itm}_\\nu^2 = 0 eV^2'),'FontWeight','normal','FontSize',get(gca,'FontSize')); 

plotdir = strrep(savedir,'results','plots');
plotname = sprintf('%sksn2_Contour_CompareRasterScanGridSearch.png',plotdir);
print(gcf,plotname,'-dpng','-r400');
fprintf('save plot to %s \n',plotname);