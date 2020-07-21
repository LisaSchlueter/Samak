T0=ref_RelicNuBkg_TDR('ToggleES','ON');
T1=ref_RelicNuBkg_TDR('ToggleES','OFF');
T0.ComputeTBDDS;
T1.ComputeTBDDS;

NormDist = 'a/sqrt(2*pi*b^2)*exp(-(x-c)^2/(2*b^2))';
pd0 = fit(T0.Te-T0.Q,T0.TBDDS_R,NormDist,'Start',[5.2e-6 0.14 -0.05]);
pd1 = fit(T1.Te-T1.Q,T1.TBDDS_R,NormDist,'Start',[5.2e-6 0.09 0]);
T0.TBDDS_R(1:200) = 0;

hT0 = semilogy((T0.Te-T0.Q),T0.TBDDS_R,'LineWidth',2,'Color','Black','LineStyle','-');
hold on;
hT1 = semilogy((T1.Te-T1.Q),T1.TBDDS_R,'LineWidth',2,'Color','Blue','LineStyle',':');
plot(pd0);
plot(pd1);
grid on;
ylim([1e-20,1e2]);
xlabel('E-E_0 (eV)','FontSize',12);
str = sprintf('dN/dE');
ylabel(str,'FontSize',14);
strT1 = sprintf('std: %1g; sigma: %1g',std(T0.TBDDS_R),pd0.b);
strT2 = sprintf('std: %1g; sigma: %1g',std(T1.TBDDS_R),pd1.b);
lh1 = legend([hT0 hT1],strT1,strT2);
%legend(lh1,'box','off');
set(lh1,'FontSize',12);
%axis([-1 2 1e-2 10000])
title(sprintf('Relic Neutrinos Capture Tritium - %.1f years',T1.TimeSec/(365.242*24*3600)),'FontSize',14);
PrettyFigureFormat;
hold off;