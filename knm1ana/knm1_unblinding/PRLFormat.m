% My pretty format for figures - G. Mention 2010.
% FontName = 'Times New Roman';
% FontSize = 18;
% AxisFontSize = 23;
FontName = 'Times New Roman';%'Helvetica';%'Times New Roman';%'Helvetica';%'Deja Vous Sans';
FontSize = 18;
AxisFontSize = 22;%24;
FontWeight = 'normal';
set(gca,'FontName',FontName,'FontSize',FontSize,'FontWeight',FontWeight)
set(get(gca,'XLabel'),'FontName',FontName,'FontSize',AxisFontSize,'FontWeight',FontWeight);
set(get(gca,'YLabel'),'FontName',FontName,'FontSize',AxisFontSize,'FontWeight',FontWeight);
set(gca,'TickLength',[.01 .01]);
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
set(gcf,'Color','w');
box on;
set(gca,'LineWidth',1.5);
set(gca,'Layer','top');
