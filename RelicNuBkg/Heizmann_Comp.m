H0 = importdata('RelicNuBkg/eta0.txt');
H1 = importdata('RelicNuBkg/eta1e11.txt');
H2 = importdata('RelicNuBkg/eta1e12.txt');
H3 = importdata('RelicNuBkg/eta1e13.txt');
I0 = importdata('RelicNuBkg/HeizmannIS_0.txt');
I1 = importdata('RelicNuBkg/HeizmannIS_1e11.txt');
I2 = importdata('RelicNuBkg/HeizmannIS_1e12.txt');
I3 = importdata('RelicNuBkg/HeizmannIS_1e13.txt');

fig1 = figure(1);
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
TBDDS_0 = T0.TBDDS.*0.008077/5e17./(0.5*(1-cos(asin(sqrt(T0.WGTS_B_T./T0.MACE_Bmax_T)))).*(T0.FPD_MeanEff*T0.FPD_Coverage) .* numel(T0.FPD_PixList)/148);
TBDDS_R_1 = T1.TBDDS_R.*0.008077/5e17./(0.5*(1-cos(asin(sqrt(T1.WGTS_B_T./T1.MACE_Bmax_T)))).*(T1.FPD_MeanEff*T1.FPD_Coverage) .* numel(T1.FPD_PixList)/148);
TBDDS_R_2 = T2.TBDDS_R.*0.008077/5e17./(0.5*(1-cos(asin(sqrt(T2.WGTS_B_T./T2.MACE_Bmax_T)))).*(T2.FPD_MeanEff*T2.FPD_Coverage) .* numel(T2.FPD_PixList)/148);
TBDDS_R_3 = T3.TBDDS_R.*0.008077/5e17./(0.5*(1-cos(asin(sqrt(T3.WGTS_B_T./T3.MACE_Bmax_T)))).*(T3.FPD_MeanEff*T3.FPD_Coverage) .* numel(T3.FPD_PixList)/148);
TBDDS_R_4 = T4.TBDDS_R.*0.008077/5e17./(0.5*(1-cos(asin(sqrt(T4.WGTS_B_T./T4.MACE_Bmax_T)))).*(T4.FPD_MeanEff*T4.FPD_Coverage) .* numel(T4.FPD_PixList)/148);

hT1 = semilogy((T1.Te-T1.Q),TBDDS_R_1,'LineWidth',2,'Color','Black','LineStyle','-');
hT2 = semilogy((T1.Te-T1.Q),TBDDS_R_2,'LineWidth',2,'Color','Red','LineStyle','-');
hT3 = semilogy((T1.Te-T1.Q),TBDDS_R_3,'LineWidth',2,'Color','Magenta','LineStyle','-');
hT4 = semilogy((T1.Te-T1.Q),TBDDS_R_4,'LineWidth',2,'Color','Green','LineStyle','-');
hT0 = semilogy((T1.Te-T1.Q),TBDDS_0,'LineWidth',2,'Color','Blue','LineStyle','-');
grid on;
ylim([1e-23,1e-18]);
xlim([-5,2]);
xlabel('E-E_0 (eV)','FontSize',12);
str = sprintf('Rate per 0.1 eV per second');
ylabel(str,'FontSize',14);
strT0 = sprintf('\\eta = 0');
strT2 = sprintf('\\eta = 1e11, R_{\\nu}^{FK} = %0.2g (R_{\\nu}^{FH} = %0.2g)',simpsons(T2.Te,TBDDS_R_2),simpsons(H1(:,1)-5,H1(:,2)));
strT3 = sprintf('\\eta = 1e12, R_{\\nu}^{FK} = %0.2g (R_{\\nu}^{FH} = %0.2g)',simpsons(T3.Te,TBDDS_R_3),simpsons(H2(:,1)-5,H2(:,2)));
strT4 = sprintf('\\eta = 1e13, R_{\\nu}^{FK} = %0.2g (R_{\\nu}^{FH} = %0.2g)',simpsons(T4.Te,TBDDS_R_4),simpsons(H3(:,1)-5,H3(:,2)));
lh1=legend([hT0 hT2 hT3 hT4],strT0,strT2,strT3,strT4);
set(lh1,'FontSize',12);
title(sprintf('Dashed lines: Florian Heizmann'))
PrettyFigureFormat;
hold off;

fig2 = figure(2);
hI0 = semilogy(I0(:,1)-5,I0(:,2),'LineWidth',2,'Color','Blue','LineStyle','--');
hold on;
hI1 = semilogy(I1(:,1)-5,I1(:,2),'LineWidth',2,'Color','Red','LineStyle','--');
hI2 = semilogy(I2(:,1)-5,I2(:,2),'LineWidth',2,'Color','Magenta','LineStyle','--');
hI3 = semilogy(I3(:,1)-5,I3(:,2),'LineWidth',2,'Color','Green','LineStyle','--');

T0.ComputeTBDIS;
T2.ComputeTBDIS;
T3.ComputeTBDIS;
T4.ComputeTBDIS;
TBDIS_0 = T0.TBDIS./T0.qUfrac./T0.TimeSec;
TBDIS_2 = T2.TBDIS./T2.qUfrac./T2.TimeSec;
TBDIS_3 = T3.TBDIS./T3.qUfrac./T3.TimeSec;
TBDIS_4 = T4.TBDIS./T4.qUfrac./T4.TimeSec;
e1 = T0.qU(T0.qU(:,1)>(T0.Q-310),:)-18575;
tmpis0 = TBDIS_0(T0.qU(:,1)>(T0.Q-310),:);
tmpis2 = TBDIS_2(T2.qU(:,1)>(T2.Q-310),:);
tmpis3 = TBDIS_3(T3.qU(:,1)>(T3.Q-310),:);
tmpis4 = TBDIS_4(T4.qU(:,1)>(T4.Q-310),:);

hTI2 = semilogy(e1,tmpis2,'-s','LineWidth',2,'Color','Red','LineStyle','-');
hTI3 = semilogy(e1,tmpis3,'-s','LineWidth',2,'Color','Magenta','LineStyle','-');
hTI4 = semilogy(e1,tmpis4,'-s','LineWidth',2,'Color','Green','LineStyle','-');
hTI0 = semilogy(e1,tmpis0,'-s','LineWidth',2,'Color','Blue','LineStyle','-');
grid on;
%ylim([1e-23,1e-18]);
xlim([-5,2]);
xlabel('E-E_0 (eV)','FontSize',12);
str = sprintf('Event Rate (Hz)');
ylabel(str,'FontSize',14);
strT0 = sprintf('\\eta = 0');
strT2 = sprintf('\\eta = 1e11, ');
strT3 = sprintf('\\eta = 1e12');
strT4 = sprintf('\\eta = 1e13');
lh1=legend([hTI0 hTI2 hTI3 hTI4],strT0,strT2,strT3,strT4);
set(lh1,'FontSize',12);
title(sprintf('Dashed lines: Florian Heizmann'))
PrettyFigureFormat;
hold off;

fig3 = figure(3);
hC1 = semilogy(e1,tmpis2-tmpis0,'LineWidth',2,'Color','Red','LineStyle','-');
hold on;
hC2 = semilogy(e1,tmpis3-tmpis0,'LineWidth',2,'Color','Magenta','LineStyle','-');
hC3 = semilogy(e1,tmpis4-tmpis0,'LineWidth',2,'Color','Green','LineStyle','-');
hF1 = semilogy(I1(:,1)-5,I1(:,2)-I0(:,2),'LineWidth',2,'Color','Red','LineStyle','--');
hF2 = semilogy(I2(:,1)-5,I2(:,2)-I0(:,2),'LineWidth',2,'Color','Magenta','LineStyle','--');
hF3 = semilogy(I3(:,1)-5,I3(:,2)-I0(:,2),'LineWidth',2,'Color','Green','LineStyle','--');
grid on;
xlim([-5,2]);
xlabel('E-E_0 (eV)','FontSize',12);
str = sprintf('Capture Event Rate (Hz)');
ylabel(str,'FontSize',14);
strT2 = sprintf('\\eta = 1e11, R_{\\nu}^{FK} = %0.2g (R_{\\nu}^{FH} = %0.2g)',tmpis2(30)-tmpis0(30),max(I1(:,2)-I0(:,2)));
strT3 = sprintf('\\eta = 1e12, R_{\\nu}^{FK} = %0.2g (R_{\\nu}^{FH} = %0.2g)',tmpis3(30)-tmpis0(30),max(I2(:,2)-I0(:,2)));
strT4 = sprintf('\\eta = 1e13, R_{\\nu}^{FK} = %0.2g (R_{\\nu}^{FH} = %0.2g)',tmpis4(30)-tmpis0(30),max(I3(:,2)-I0(:,2)));
lh1=legend([hC1 hC2 hC3],strT2,strT3,strT4);
set(lh1,'FontSize',12);
title(sprintf('Dashed lines: Florian Heizmann'))
PrettyFigureFormat;
hold off;