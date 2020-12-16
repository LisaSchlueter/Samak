
load('CATSfit_Hat_Qsq_Uniform_CMShape.mat');

GetFigure;

p4 = plot(qU(11:end)-18574,Hat(:,4),':','LineWidth',2,'Color',rgb('LimeGreen'));
hold on;
p3 = plot(qU(11:end)-18574,Hat(:,3),'--','LineWidth',2,'Color',rgb('DodgerBlue'));
p1 = plot(qU(11:end)-18574,Hat(:,1),'-','LineWidth',2.5,'Color',rgb('Orange'));
hold on;
p2 = plot(qU(11:end)-18574,Hat(:,2),'-.','LineWidth',2,'Color',rgb('Red'));

p5 = plot(qU(11:end)-18574,sum(Hat,2),'LineWidth',2.5,'Color',rgb('Gray'));
leg =legend([p1,p2,p3,p4,p5],sprintf('{\\itm}_\\nu^2'),sprintf('{\\itE}_0^{fit}'),...
    sprintf('{\\itB}'),sprintf('{\\itN}'),'Total');%,'{\\itB}','{\\itN}',
PrettyLegendFormat(leg);
PrettyFigureFormat;
ylabel('Leverage')
xlabel('Retarding potential - 18574 (eV)');

print(gcf,'plots/LeveragemNuE0CM.png','-dpng','-r350');