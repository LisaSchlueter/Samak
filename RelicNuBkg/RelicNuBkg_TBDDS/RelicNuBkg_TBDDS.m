
ToggleES = 'ON';
mnuSq_i = 0;
eta_i = 1;
Time = 1000*24*3600;

% init TBD object with Design Report (TDR) properties:
T1 = ref_RelicNuBkg_TDR('ToggleES',ToggleES,'mnuSq_i',mnuSq_i,'eta_i',eta_i,'TimeSec',Time);
T2 = ref_RelicNuBkg_TDR('ToggleRelic','OFF','mnuSq_i',mnuSq_i,'TimeSec',Time);

% Compute diff spec.
T1.ComputeTBDDS;
T2.ComputeTBDDS;
T1.TBDDS = T1.TBDDS;
T2.TBDDS = T2.TBDDS;

hT1 = semilogy((T1.Te-T1.Q),T1.TBDDS,'LineWidth',2,'Color','Black','LineStyle','-');
hold on;
hT2 = semilogy((T1.Te-T1.Q),T2.TBDDS,'LineWidth',2,'Color','Blue','LineStyle',':');
hA1 = semilogy((T1.Te-T1.Q),T1.TBDDS_R,'LineWidth',2,'Color','Red','LineStyle','--');
grid on;
ylim([1e-20,1e2]);
xlabel('E-E_0 (eV)','FontSize',12);
str = sprintf('dN/dE');
ylabel(str,'FontSize',14);
strT1 = sprintf('Tritium Beta Decay Spectrum');
strT2 = sprintf('Beta Decay + E-capture');
str1 = sprintf('E-capture: %.3g evts',T1.NormFactorTBDDS_R*0.57);
lh1 = legend([hT1 hT2 hA1],strT1,strT2,str1);
%legend(lh1,'box','off');
set(lh1,'FontSize',12);
%axis([-1 2 1e-2 10000])
title(sprintf('Relic Neutrinos Capture Tritium - %.1f years',T1.TimeSec/(365.242*24*3600)),'FontSize',14);
PrettyFigureFormat;
hold off;