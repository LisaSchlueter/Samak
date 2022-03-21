% My pretty format for figures - G. Mention 2010.
function PrettyFigureFormat(varargin)
p=inputParser;
p.addParameter('FontSize',20,@(x)isfloat(x));
p.parse(varargin{:});
FontSize = p.Results.FontSize;

FontName = 'Arial';%'Helvetica';%'Helvetica';%'Helvetica';%'Deja Vous Sans';
FontWeight = 'normal';
set(gca,'FontName',FontName,'FontSize',FontSize-4,'FontWeight',FontWeight)
set(get(gca,'XLabel'),'FontName',FontName,'FontSize',FontSize,'FontWeight',FontWeight);
set(get(gca,'YLabel'),'FontName',FontName,'FontSize',FontSize,'FontWeight',FontWeight);
set(gca,'TickLength',[.01 .01]);
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
set(gcf,'Color','w');
box on;
set(gca,'LineWidth',1.5);
set(gca,'Layer','top');
end