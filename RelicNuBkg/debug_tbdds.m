T = ref_RelicNuBkg_TDR;
T.ComputeTBDDS;
hT = semilogy((T.Te-T.Q),T.TBDDS,'LineWidth',2,'Color','Blue','LineStyle',':');
grid on;
xlabel('E-E_0 (eV)','FontSize',12);
str = sprintf('Rate per energy bin');
ylabel(str,'FontSize',14);
strT = sprintf('Tritium Beta Decay Spectrum');
lh1 = legend([hT],str);
%legend(lh1,'box','off');
set(lh1,'FontSize',12);
PrettyFigureFormat;