H0 = importdata('RelicNuBkg/FormaggioISBeta.txt');
H1 = importdata('RelicNuBkg/FormaggioISRelic.txt');

T = RelicNuDebug('Params','Formaggio','ToggleES','ON');
T.SetRhoDmnu2TimeEta('mnu',1,'Time',365.242*24*3600,'Eta',1);%,'RhoDMode','EffMass','RhoD',66.5e-6);
fig1=figure(1);
T.TotalCountsStudy('Plot','ON');

%U = RelicNuDebug('Params','TDR','ToggleES','OFF');
%U.SetRhoDmnu2TimeEta('mnu',1,'Time',365.242*24*3600,'Eta',1,'RhoDMode','EffMass','RhoD',66.5e-6);
%U.SetMass('EffMass','ON','tmass',66.5e-6);
%U.TotalCountsStudy;

%Eta = 1e6./(T.R.NormFactorTBDDS_R*T.R.WGTS_CD_MolPerCm2/5e17./T.R.TimeSec.*365.242.*24.*3600);
%T.SetRhoDmnu2TimeEta('mnu',1,'Time',365.242*24*3600,'Eta',Eta);

T.SetRhoDmnu2TimeEta('mnu',1,'Time',365.242*24*3600,'Eta',(6.4915/5)*1.7931e11/0.5717);%,'RhoDMode','EffMass','RhoD',66.5e-6);
T.TotalCountsStudy('Plot','OFF');
fig2=figure(2);
hH0 = semilogy(H0(:,1)-T.R.Q,H0(:,2),'-s','LineWidth',2,'Color','Green','LineStyle','--');
hold on;
hH1 = semilogy(H1(:,1)-T.R.Q,H1(:,2),'-s','LineWidth',2,'Color','Red','LineStyle','--');
TBDIS_R = T.R.TBDIS./T.R.qUfrac./T.R.TimeSec;
TBDIS_B = (T.R.TBDIS-T.R.TBDIS_R)./T.R.qUfrac./T.R.TimeSec;
           
e1 = T.R.qU(T.R.qU(:,1)>(T.R.Q-310),:)-18575;
tmpis1 = TBDIS_B(T.R.qU(:,1)>(T.R.Q-310),:);
tmpis2 = TBDIS_R(T.R.qU(:,1)>(T.R.Q-310),:);

hT1 = semilogy((e1),tmpis1,'-s','LineWidth',2,'Color','Black','LineStyle','-');
hT2 = semilogy((e1),tmpis2,'-s','LineWidth',2,'Color','Blue','LineStyle','-');

grid on;
%ylim([1e-20,1e2]);
xlabel('E-E_0 (eV)','FontSize',12);
str = sprintf('Event Rate (Hz)');
ylabel(str,'FontSize',14);
strT1 = sprintf('0 capture events/yr');
strT2 = sprintf('1e6 capture events/yr');
strH0 = sprintf('Formaggio 0 capture events/yr');
strH1 = sprintf('Formaggio 1e6 capture events/yr');
lh1 = legend([hT1 hT2 hH0 hH1],strT1,strT2,strH0,strH1);
set(lh1,'FontSize',12);
PrettyFigureFormat;
hold off;

fig3=figure(3);
p=plot(e1,(H0(:,2)-tmpis1)./tmpis1,'-s','LineWidth',2,'Color','Black','LineStyle','-');
grid on;
xlabel('E-E_0 (eV)','FontSize',12);
ylabel('Residual','FontSize',14);
str = sprintf('Residuals for 0 capture events, column density 5e17');
lh1 = legend([p], str);
set(lh1,'FontSize',12);
PrettyFigureFormat;


T.SetRhoDmnu2TimeEta('mnu',1,'Time',365.242*24*3600,'Eta',1.7931e11/0.5717,'RhoDMode','EffMass','RhoD',66.5e-6);
T.TotalCountsStudy('Plot','OFF');
fig4=figure(4);
hH0 = semilogy(H0(:,1)-T.R.Q,H0(:,2),'-s','LineWidth',2,'Color','Green','LineStyle','--');
hold on;
hH1 = semilogy(H1(:,1)-T.R.Q,H1(:,2),'-s','LineWidth',2,'Color','Red','LineStyle','--');
TBDIS_R = T.R.TBDIS./T.R.qUfrac./T.R.TimeSec;
TBDIS_B = (T.R.TBDIS-T.R.TBDIS_R)./T.R.qUfrac./T.R.TimeSec;
           
e1 = T.R.qU(T.R.qU(:,1)>(T.R.Q-310),:)-18575;
tmpis1 = TBDIS_B(T.R.qU(:,1)>(T.R.Q-310),:);
tmpis2 = TBDIS_R(T.R.qU(:,1)>(T.R.Q-310),:);

hT1 = semilogy((e1),tmpis1,'-s','LineWidth',2,'Color','Black','LineStyle','-');
hT2 = semilogy((e1),tmpis2,'-s','LineWidth',2,'Color','Blue','LineStyle','-');

grid on;
%ylim([1e-20,1e2]);
xlabel('E-E_0 (eV)','FontSize',12);
str = sprintf('Event Rate (Hz)');
ylabel(str,'FontSize',14);
strT1 = sprintf('0 capture events/yr');
strT2 = sprintf('1e6 capture events/yr');
strH0 = sprintf('Formaggio 0 capture events/yr');
strH1 = sprintf('Formaggio 1e6 capture events/yr');
lh1 = legend([hT1 hT2 hH0 hH1],strT1,strT2,strH0,strH1);
set(lh1,'FontSize',12);
PrettyFigureFormat;
hold off;

fig5=figure(5);
p=plot(e1,(H0(:,2)-tmpis1)./tmpis1,'-s','LineWidth',2,'Color','Black','LineStyle','-');
grid on;
xlabel('E-E_0 (eV)','FontSize',12);
ylabel('Residual','FontSize',14);
str = sprintf('Residuals for 0 capture events, large column density');
lh1 = legend([p], str);
set(lh1,'FontSize',12);
PrettyFigureFormat;