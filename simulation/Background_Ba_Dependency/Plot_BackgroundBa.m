% Compute Background
Ba = 1e-04.*(3:0.1:12);
WGTS_B_T = 3.6*0.7; 
BKG70_Nikolaus  = 1e3.*GetBackground('MACE_Ba_T',Ba,'WGTS_B_T',WGTS_B_T,'AnchorBkg6G',0.2376);
BKG70_FT = 1e3.*GetBackground('MACE_Ba_T',Ba,'WGTS_B_T',WGTS_B_T,'FPD_ROIlow',14);
BKG70_FT_ROI26 = 1e3.*GetBackground('MACE_Ba_T',Ba,'WGTS_B_T',WGTS_B_T,'FPD_ROIlow',26);


%% plot background and resolution
close;
f24 = figure('Renderer','openGL');
set(f24, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.8]);
PrettyFigureFormat;
Bmax = 6*0.7; % pinch magnet
R70 = 18575*Ba/Bmax;

[ax, h1, h2]  = plotyy(Ba*1e4,BKG70_FT,Ba*1e4,R70); hold on;%BKG70_Nikolaus
%set(ax,'NextPlot','add');
%p70FT = plot(ax(1),Ba*1e4,BKG70_FT,'-.','Color',rgb('CadetBlue'),'LineWidth',4);
%p70FT_26ROI = plot(ax(1),Ba*1e4,BKG70_FT_ROI26,':','Color',rgb('CadetBlue'),'LineWidth',4);
set(h1,'LineWidth',4,'Color',rgb('CadetBlue')); set(h2,'LineWidth',4,'Color',rgb('IndianRed'));%,'LineStyle','none');
set(ax,'FontWeight','Bold','LineWidth',2);
set(ax,{'ycolor'},{rgb('CadetBlue'); rgb('IndianRed')});%
set(get(ax(1),'Ylabel'),'String','background rate 148 pixels (mcps)','FontWeight','Bold');
set(get(ax(2),'Ylabel'),'String','resolution \DeltaE (eV)','FontWeight','Bold');
set(ax,'FontSize',20)
xlabel('magnetic field analyzing plane B_a (10^{-4} T)');
set(ax,'XLim',1e4.*[min(Ba) max(Ba)]);
set(ax(1),'YTick',[0:100:950]);
set(ax(2),'YTick',[0:1:6]);
set(ax(1),'YLim',[0 950]);%[0,750]);
set(ax(2),'YLim',[0,6]);
%set(ax(2),'YTickLabel',''); set(ax(2),'YTick',[]);
grid on;
%leg = legend([h1,p70FT_26ROI,p70FT],'SDS measurements','First Tritium Campaign (26 keV ROI cut)','First Tritium Campaign (14 keV ROI cut)');
%%leg.Title.String = 'Background Measurements';
leg = legend(h1,'First Trititum Campaign (14 keV ROI cut)');
leg.Location = 'southwest';
legend boxoff;
%title(sprintf('Nominal KATRIN settings \nB_s = 2.5T  -  B_{max} = 4.2T (70%%)'));
save_name = '../../simulation/Background_Ba_Dependency/plots/BackgroundBaResolution';
print([save_name,'.png'],'-dpng','-r350');
%publish_figurePDF(f24,[save_name,'.pdf']);

 %% plot background only
% f23 = figure(23);
% set(f23, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.7]);
% p70 = plot(Ba*1e4,BKG70_Nikolaus,'-','LineWidth',3,'Color',rgb('CadetBlue'));
% xlabel('B_a (10^{-4} T)');
% ylabel('Background Rate (mcps)');
% xticks([3:1:12]);
% xlim([3 12]);
% ylim([0 750]);
% set(gca,'XScale','log');
% PrettyFigureFormat;
% set(gca,'FontSize',16)
% grid on;
% hold on;
% %p100 = plot(Ba*1e4,BKG100*1e3,'--o','LineWidth',3,'Color',rgb('IndianRed'));
% %legend([p70,p100],'Bs = 2.54 T', 'Bs = 3.60 T');
% %print(gcf,'../../simulation/Background_Ba_Dependency/plots/Background_Ba_dependence.png','-dpng');
% title(sprintf('Nominal KATRIN settings \nB_s = 2.5T  -  B_{max} = 4.2T (70%%)'));
% %print('../../simulation/Background_Ba_Dependency/plots/BackgroundBa_70percentBFields.png','-dpng');
