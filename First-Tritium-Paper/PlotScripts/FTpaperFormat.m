% My pretty format for figures - G. Mention 2010.

FontName = 'Arial';
FontSize = 19;
FontWeight = 'normal';
set(gca,'FontName',FontName,'FontSize',FontSize,'FontWeight',FontWeight)
set(get(gca,'XLabel'),'FontName',FontName,'FontSize',FontSize,'FontWeight',FontWeight);
set(get(gca,'YLabel'),'FontName',FontName,'FontSize',FontSize,'FontWeight',FontWeight);
set(gca,'TickLength',[.01 .01]);
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
set(gcf,'Color','w');
box on;
orange=1/255*[255,140,0];
myblue = rgb('CadetBlue');
%set(gca,'LineWidt',2.5);
set(gca,'LineWidt',1);
set(gca,'Layer','top')