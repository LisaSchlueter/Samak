A = ref_TBD_NominalKATRIN;
A.ComputeTBDDS;

%%
close;
f11 = figure('Renderer','openGL');
set(f11, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.7]);
plot(A.TTexE,A.TTexP*100,'LineWidth',5,'color',rgb('CadetBlue'));
set(gca,'XScale','log');
xlabel('excitation energy (eV)');
ylabel('probability (%)');
xlim([0.1 200]);
leg = legend('T_2 (Saenz et. al.)'); legend boxoff
leg.FontSize = 18;
PrettyFigureFormat;
set(gca,'FontSize',18);
save_name = 'TT-FSD_TheorySaenz';
print(f11,['./plots/',save_name,'.png'],'-dpng');