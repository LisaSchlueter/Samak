% test script
% can be deleted later

% ini TBD object with Design Report (TDR) properties:
A = ref_RelicNuBkg_TDR('DopplerEffectFlag','FSD');
A.DE_sigma = 0.0935;
A.LoadFSD;
% Compute diff spec.
A.ComputeTBDDS;
A.TBDDS = A.TBDDS./simpsons(A.Te,A.TBDDS);
pFSD = A.PlotTBDDS;
hold on;

A.TTFSD = 'OFF';
A.DTFSD = 'OFF';
A.HTFSD = 'OFF';

A.ComputeTBDDS;
pNoFSD = plot(A.Te-A.Q,A.TBDDS./simpsons(A.Te,A.TBDDS),'-.','LineWidth',3,'Color',rgb('Orange'));
hold off;
%%
leg = legend([pNoFSD,pFSD],'Without FSD','With FSD');
leg.Title.String = sprintf('Differential \\beta-spectrum');
leg.Title.FontWeight = 'normal';
xlim([-15 1]);
set(gca,'YScale','log');
ylim([5e-100 10]);
