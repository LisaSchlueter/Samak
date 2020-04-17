
% bar plot for background systematic
PlotColor = {rgb('DodgerBlue'),rgb('PowderBlue'),rgb('FireBrick'),rgb('DarkOrange')};
SingleBarY = [0.085,0.048,0.080,0.030];

bsingle = cell(numel(SingleBarY),1);
t =cell(numel(SingleBarY)+1,1);
SingleBarX = numel(SingleBarY):-1:1;
figure('Units','normalized','Position',[0.1,0.1,0.65,0.6]);
for i=1:numel(SingleBarY)
    bsingle{i}  = barh(SingleBarX(i),SingleBarY(i));
    hold on;
    bsingle{i}.FaceColor = PlotColor{i}; bsingle{i}.LineStyle ='none';
    bsingle{i}.FaceAlpha=1;
    bsingle{i}.BarWidth = 0.9;
    
    if SingleBarY(i)<0.001
        tstr = sprintf('<10^{-3}') ;
    elseif SingleBarY(i)<0.01 && SingleBarY(i)>0.001
        tstr = sprintf('%.0f\\cdot10^{-3}',SingleBarY(i)*1e3) ;
    else
        tstr = sprintf('%.1f\\cdot10^{-2}',SingleBarY(i)*1e2) ;
    end
    t{i}= text(0.098,SingleBarX(i),tstr,...%max(SingleBarY)*1.4
        'HorizontalAlignment','right','FontSize',20,...
        'Interpreter','tex');
end

xlabel(sprintf('Sensitivity {\\itm}^2_\\nu at %.0f%% C.L. (eV^{ 2})',0.683*100));
xmin = 0;
xmax = 0.1;% max(SingleBarY)*2;
xlim([xmin, xmax]);

PrettyFigureFormat('FontSize',22);
% no y-ticks
set(gca,'YMinorTick','off');
set(gca,'XMinorTick','off');
h = gca;
h.YRuler.TickLength = [0,0];
yticks(flip(SingleBarX))
yticklabels({'Slope fit + FT','Slope fit unconstrained','Gauss + FT ','Gauss + KNM1/2'});
title('Background slope systematic uncertainty','FontWeight','normal');

savedir = [getenv('SamakPath'),'knm2ana/knm2_systematics/plots/'];
MakeDir(savedir);
savename = sprintf('%sKNM2BkgSysBreakdown.png',savedir);
print(gcf,savename,'-dpng','-r400');