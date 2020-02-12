
% FPD Coverage versus ROI
% Sanshiro data, 31/10/2018
fpdF = load('FPDROICoverage.txt');
fpdCov = @(ROIdown) interp1(fpdF(:,1),fpdF(:,2),ROIdown);

% Background coverage versus ROI
% Anna data, 31/10/2018
bkgF = load('BackgroundROICut.txt');
bkgCov = @(ROIdown) interp1(bkgF(:,1),bkgF(:,2)*1000,ROIdown);

% Plot
ROIlow = 14:.3:26;

clf
figure(1)
set(1, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.7]);
yyaxis left

title('[FPD Coverage , Background] Versus ROI')

hold on

xlabel('FPD ROI Lower Limit (keV)')
ylabel('FPD Coverage (FT wafer)')

yyaxis right
ylabel('Background 148 pixels (mcps)')

yyaxis left
hl=plot(ROIlow,smooth(fpdCov(ROIlow)),'LineWidth',8);

yyaxis right
hr=plot(ROIlow,smooth(bkgCov(ROIlow)),'LineWidth',8);

hold off

grid on

FontName = 'Helvetica';
FontSize = 28;
FontWeight = 'bold';
set(gca,'FontName',FontName,'FontSize',FontSize,'FontWeight',FontWeight)
set(get(gca,'XLabel'),'FontName',FontName,'FontSize',FontSize,'FontWeight',FontWeight);
set(get(gca,'YLabel'),'FontName',FontName,'FontSize',FontSize,'FontWeight',FontWeight);
set(gca,'TickLength',[.02 .02]);
set(gca,'XMinorTick','off');
set(gca,'YMinorTick','off');
set(gcf,'Color','w');
box on;
orange=1/255*[255,140,0];
myblue = rgb('CadetBlue');

%%
figure(2)
set(2, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.9, 0.7]);
%subplot(3,1,[1 2])
yyaxis left

%title('[FPD Coverage , Background] Versus ROI')

hold on

xlabel('FPD ROI Lower Limit (keV)')
ylabel('Tritium Signal Loss (%)')

yyaxis right
ylabel('Background 148 pixels (mcps)')

yyaxis left
hl=plot(ROIlow,(1-smooth(fpdCov(ROIlow)))*100,'LineWidth',8);

yyaxis right
hr=plot(ROIlow,smooth(bkgCov(ROIlow)),'LineWidth',8);

hold off

grid on

FontName = 'Helvetica';
FontSize = 34;
FontWeight = 'bold';
set(gca,'FontName',FontName,'FontSize',FontSize,'FontWeight',FontWeight)
set(get(gca,'XLabel'),'FontName',FontName,'FontSize',FontSize,'FontWeight',FontWeight);
set(get(gca,'YLabel'),'FontName',FontName,'FontSize',FontSize,'FontWeight',FontWeight);
set(gca,'TickLength',[.02 .02]);
set(gca,'XMinorTick','off');
set(gca,'YMinorTick','off');
set(gcf,'Color','w');
box on;
orange=1/255*[255,140,0];
myblue = rgb('CadetBlue');
publish_figurePDF(2,'FPDCoverageBackgroundVersusROI.pdf')
% subplot(3,1,3)
% hold on
% h2l=plot(ROIlow,smooth(1000*fpdCov(ROIlow)./sqrt(bkgCov(ROIlow))),'LineWidth',5,'Color',rgb('Black'),'LineStyle','--');
% ylabel('10cps /sqrt(B)')
% hold off
% 
% FontName = 'Helvetica';
% FontSize = 28;
% FontWeight = 'bold';
% set(gca,'FontName',FontName,'FontSize',FontSize,'FontWeight',FontWeight)
% set(get(gca,'XLabel'),'FontName',FontName,'FontSize',FontSize,'FontWeight',FontWeight);
% set(get(gca,'YLabel'),'FontName',FontName,'FontSize',FontSize,'FontWeight',FontWeight);
% set(gca,'TickLength',[.02 .02]);
% set(gca,'XMinorTick','off');
% set(gca,'YMinorTick','off');
% set(gcf,'Color','w');
% box on;
% orange=1/255*[255,140,0];
% myblue = rgb('CadetBlue');