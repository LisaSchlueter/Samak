ELossFlag = 'Abdurashitov';
A = ref_TBD_NominalKATRIN('TD','Flat100','ELossFlag',ELossFlag);
f = A.ComputeELossFunction;
%%
maxE = (A.Q_i-A.qUmin)*1.5; minE=-maxE; NbinE = (maxE-minE)/0.5;
E = linspace(minE,maxE,NbinE);

f1 = f{1};
f2 = f{2};
f3 = f{3};
f4 = f{4};
%f5 = f{5};
%f6 = f{6};


fig5 = figure('Renderer','opengl');
set(fig5, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7 ,0.7]);
plot(E,f1(E),'Linewidth',3,'Color',rgb('FireBrick'));
hold on;
plot(E,f2(E),'Linewidth',3,'Color',rgb('Orange'),'LineStyle',':');
plot(E,f3(E),'Linewidth',3,'Color',rgb('RoyalBlue'),'LineStyle','--');
plot(E,f4(E),'Linewidth',3,'Color',rgb('CadetBlue'),'LineStyle','-.');
%plot(E,f5(E),'Linewidth',3,'Color',rgb('Black'),'LineStyle','-.');
%plot(E,f6(E),'Linewidth',3,'Color',rgb('CadetBlue'),'LineStyle','-.');

PrettyFigureFormat;
set(gca,'FontSize',24);
xlim([0 80])
ylim([0 0.27])
xlabel(['energy loss ',char(949),' (eV)']);
ylabel(['probability f(',char(949),')']);
leg = legend(string(1:4)); legend boxoff;
leg.FontSize = 24;
leg.Title.String = 'number of scatterings';

save_name = sprintf('../../studies/ft_CovarianceMatrices/plots/EnergyLossFunctions/ELossFunctions_%s',ELossFlag);
print(fig5,[save_name,'.png'],'-dpng','-r300');
publish_figurePDF(fig5,[save_name,'.pdf']);

