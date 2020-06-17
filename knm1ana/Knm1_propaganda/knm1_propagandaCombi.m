

fcombi = figure('Units','normalized','Position',[0.1,0.1,0.5,1]);
s1 = subplot(3,1,1);
[leg1,ax1,zoomBox] = m2resultsVsyear2019('Combi','ON','LegPos','northoutside','SavePlot','OFF','LegOpt','Short');

mypos1 = ax1.Position;
%ax1.Position = [mypos1(1) mypos1(2)-0.12 mypos1(3) mypos1(4)+0.15];
ax1.XTickLabel = '';
ax1.XLabel.String = '';%(' ');

%%
s2 = subplot(3,1,2);
m2sysVsyear2019('Combi','ON','SavePlot','OFF')
ax2 = gca;
mypos2 = ax2.Position;
xticklabels('');
xlabel('');
ylabel({sprintf('{\\itm}_\\nu^2 systematic');sprintf(' uncertainty (eV^2)')});
set(gca,'TickDir','out');

grid off;
%% subplot 3
s3 = subplot(3,1,3);
m2statVsyear2019('Combi','ON','SavePlot','OFF')
ax3 = gca;
mypos3 = ax3.Position;
grid off;
ylabel({sprintf('{\\itm}_\\nu^2 statistical');sprintf(' uncertainty (eV^2)')});
set(gca,'TickDir','out');
%% shift plots: make space for legend and reduce white space
ax1.Position = [mypos1(1)+0.05 mypos1(2)-0.125 mypos1(3) mypos1(4)+0.2];
ax2.Position = [mypos2(1)+0.05 mypos2(2)-0.072 mypos2(3) mypos2(4)+0.023];
ax3.Position = [mypos3(1)+0.05 mypos3(2)-0.02 mypos3(3) mypos3(4)+0.023];

zoomBox.Position = [0.542 0.605 .393 .16];% [0.542 0.58 .393 .16];

zoomBox.Legend.delete
zoomBox.YLim= [-4.5 3.5];
zoomBox.YTick = [-2 0 2];
zoomBox.XTick = [2005 2011 2019];
zoomBox.FontSize = 22;
ax1.YLabel.FontSize = 24.5;
ax2.YLabel.FontSize = 24.5;
ax3.YLabel.FontSize = 24.5;
ax3.XLabel.FontSize = 24.5;

zoomBox.XAxis.FontSize = ax1.XAxis.FontSize;
zoomBox.YAxis.FontSize = ax1.YAxis.FontSize;
zoomBox.YLabel.FontSize = ax1.YLabel.FontSize;
ax2.XTick = [];
ax2.YTick= [ 1 10 1e2];
ax3.YTick= [ 1 10 1e2];
ax2.YLabel.Position(1) = 1.9867e+03;
ax2.YLabel.Position(2) = 5.5;
ax3.YLabel.Position(1) = 1.9867e+03;
ax1.YLabel.Position(1) = 1.986e+03;
linkaxes([s1,s2,s3],'x')
ax1.XLim = [1.98951e3 2021];
ax1.YLim = [-240 70];
ax2.YLim = [0.1 500];
ax3.YLim = [0.1 500];
zoomBox.XLim = [2004.2 2020.2];
set(gca,'TickDir','out');

leg1.FontSize = 15.4;
leg1.Position = [0.122   0.91    0.88    0.0916];%.895
%% text abc)
text(1990,1.9e11,'a)','FontSize',get(gca,'FontSize')+4,'FontName',get(gca,'FontName'));   
text(1990,1.2e6,'b)','FontSize',get(gca,'FontSize')+4,'FontName',get(gca,'FontName'));            
text(1990,230,'c)','FontSize',get(gca,'FontSize')+4,'FontName',get(gca,'FontName'));            

%%
savedir = [getenv('SamakPath'),'knm1ana/Knm1_propaganda/plots/'];
MakeDir(savedir);
savename = sprintf('%sknm1_m2Vsyear2019Combi',savedir);
print(gcf,[savename,'.png'],'-dpng','-r300');
export_fig(gcf,[savename,'.pdf']);