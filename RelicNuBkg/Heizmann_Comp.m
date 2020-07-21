H0 = importdata('RelicNuBkg/eta0.txt');
H1 = importdata('RelicNuBkg/eta1e11.txt');
H2 = importdata('RelicNuBkg/eta1e12.txt');
H3 = importdata('RelicNuBkg/eta1e13.txt');

hH0 = semilogy(H0(:,1)-5,H0(:,2),'LineWidth',2,'Color','Blue','LineStyle','--');
hold on;
hH1 = semilogy(H1(:,1)-5,H1(:,2),'LineWidth',2,'Color','Red','LineStyle','--');
hH2 = semilogy(H2(:,1)-5,H2(:,2),'LineWidth',2,'Color','Magenta','LineStyle','--');
hH3 = semilogy(H3(:,1)-5,H3(:,2),'LineWidth',2,'Color','Green','LineStyle','--');

CD = 5e17;%1./(2*pi*4.5^2*0.973); %column density: one tritium atom

T0=ref_RelicNuBkg_DesignReport('mnuSq_i',1,'eta_i',0,'TimeSec',1,'WGTS_CD_MolPerCm2',CD);
T1=ref_RelicNuBkg_DesignReport('mnuSq_i',1,'eta_i',1,'TimeSec',1,'WGTS_CD_MolPerCm2',CD,'ToggleES','OFF');
T2=ref_RelicNuBkg_DesignReport('mnuSq_i',1,'eta_i',1e11,'TimeSec',1,'WGTS_CD_MolPerCm2',CD,'ToggleES','OFF');
T3=ref_RelicNuBkg_DesignReport('mnuSq_i',1,'eta_i',1e12,'TimeSec',1,'WGTS_CD_MolPerCm2',CD,'ToggleES','OFF');
T4=ref_RelicNuBkg_DesignReport('mnuSq_i',1,'eta_i',1e13,'TimeSec',1,'WGTS_CD_MolPerCm2',CD,'ToggleES','OFF');
T0.ComputeTBDDS;
T1.ComputeTBDDS;
T2.ComputeTBDDS;
T3.ComputeTBDDS;
T4.ComputeTBDDS;
T0.TBDDS = T0.TBDDS.*0.008077/5e17;
T1.TBDDS_R = T1.TBDDS_R.*0.008077/5e17;
T2.TBDDS_R = T2.TBDDS_R.*0.008077/5e17;
T3.TBDDS_R = T3.TBDDS_R.*0.008077/5e17;
T4.TBDDS_R = T4.TBDDS_R.*0.008077/5e17;
%T0.TBDDS_R = T0.TBDDS_R.*numel(T0.TBDDS)./30; %per eV
%T1.TBDDS_R = T1.TBDDS_R.*numel(T1.TBDDS)./30;
%T2.TBDDS_R = T2.TBDDS_R.*numel(T2.TBDDS)./30;
%T3.TBDDS_R = T3.TBDDS_R.*numel(T3.TBDDS)./30;
%T4.TBDDS_R = T4.TBDDS_R.*numel(T4.TBDDS)./30;
%T1.TBDDS_R = T1.TBDDS_R.*0.1;
%T2.TBDDS_R = T2.TBDDS_R.*0.1;
%T3.TBDDS_R = T3.TBDDS_R.*0.1;
%T4.TBDDS_R = T4.TBDDS_R.*0.1;
hT1 = semilogy((T1.Te-T1.Q),T1.TBDDS_R,'LineWidth',2,'Color','Black','LineStyle','-');
hT2 = semilogy((T1.Te-T1.Q),T2.TBDDS_R,'LineWidth',2,'Color','Red','LineStyle','-');
hT3 = semilogy((T1.Te-T1.Q),T3.TBDDS_R,'LineWidth',2,'Color','Magenta','LineStyle','-');
hT4 = semilogy((T1.Te-T1.Q),T4.TBDDS_R,'LineWidth',2,'Color','Green','LineStyle','-');
hT0 = semilogy((T1.Te-T1.Q),T0.TBDDS,'LineWidth',2,'Color','Blue','LineStyle','-');
grid on;
ylim([1e-23,1e-18]);
xlim([-5,2]);
xlabel('E-E_0 (eV)','FontSize',12);
str = sprintf('Rate per energy bin');
ylabel(str,'FontSize',14);
strT0 = sprintf('eta = 0');
strT2 = sprintf('eta = 1e11, events: %0.2g (%0.2g)',simpsons(T2.Te,T2.TBDDS_R),simpsons(H1(:,1)-5,H1(:,2)));
strT3 = sprintf('eta = 1e12, events: %0.2g (%0.2g)',simpsons(T3.Te,T3.TBDDS_R),simpsons(H2(:,1)-5,H2(:,2)));
strT4 = sprintf('eta = 1e13, events: %0.2g (%0.2g)',simpsons(T4.Te,T4.TBDDS_R),simpsons(H3(:,1)-5,H3(:,2)));
lh1=legend([hT0 hT2 hT3 hT4],strT0,strT2,strT3,strT4);
set(lh1,'FontSize',12);
title(sprintf('Dashed lines: Florian Heizmann'))
PrettyFigureFormat;
hold off;