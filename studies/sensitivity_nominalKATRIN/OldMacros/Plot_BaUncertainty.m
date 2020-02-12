
Ba     = [3:12]*1e-04;
BaErr = (Ba.*3.9*1e-03+1.07*1e-06);
S = GetSysBudget('SysBudget','06');

f24 = figure('Renderer','openGL');
set(f24, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.75, 0.8]);

pnone = plot(3:12,NaN.*zeros(10,1),'Color',rgb('White')); hold on;
[ax, h1, h2]  = plotyy(Ba.*1e4,BaErr*1e4,Ba.*1e4,(BaErr./Ba)*100); 
%PrettyFigureFormat;
set(ax,'XLim',[3,12]);
set(ax,'FontWeight','Bold','LineWidth',3);
set(ax,{'ycolor'},{rgb('FireBrick');rgb('CadetBlue')});
set(h1,'LineWidth',5,'Color',rgb('FireBrick'));
set(h2,'LineWidth',5,'Color',rgb('CadetBlue'));
set(get(ax(1),'Ylabel'),'String','absolute B_a uncertainty (G)','FontWeight','Bold');
set(get(ax(2),'Ylabel'),'String','relative B_a uncertainty (%)','FontWeight','Bold');
keg = legend(pnone,['\Delta B_a=',S.MACE_Ba_T_Err_str],'Location','north');
legend boxoff
leg.FontSize = 22;
set(ax,'FontSize',22);
xlim([3,12]);
xticks(3:12)
xlabel('B_a (G)')
set(ax,'TickLength',[.01 .01]);
set(ax,'XMinorTick','on');
set(ax,'YMinorTick','on');
set(gcf,'Color','w');
box on;

save_name = '../ft_CovarianceMatrices/plots/Ba_Uncertainty';
print([save_name,'.png'],'-dpng');
publish_figurePDF(f24,[save_name,'.pdf']);