T = ref_RelicNuBkg_TDR;
T.ComputeTBDDS('NormFlag','ON');
%% plot
GetFigure;
hT = semilogy((T.Te-T.Q),T.TBDDS,'LineWidth',2,'Color','Blue','LineStyle','-.');
hold on;
plot((T.Te(T.TBDDS-T.Q>0)-T.Q),T.TBDDS(T.TBDDS-T.Q>0),'LineWidth',2,'Color',rgb('Orange'))
xlabel('E-E_0 (eV)','FontSize',12);
str = sprintf('Rate per energy bin');
ylabel(str,'FontSize',14);
strT = sprintf('Tritium Beta Decay Spectrum');
lh1 = legend([hT],str);
%legend(lh1,'box','off');
lh1.FontSize = 20;
lh1.Location = 'northwest';
lh1.EdgeColor = rgb('Silver');
PrettyFigureFormat;
ylim([1e-40,20])