H0 = importdata('RelicNuBkg/FormaggioISBeta.txt');
H1 = importdata('RelicNuBkg/FormaggioISRelic.txt');

T = RelicNuAnalysis('Params','TDR','ToggleES','ON');
T.SetRhoDmnu2TimeEta('mnu',1,'Time',365.242*24*3600,'Eta',1);%,'RhoDMode','EffMass','RhoD',66.5e-6);
fig1=figure(1);
T.TotalCountsStudy('Plot','ON');

U = RelicNuAnalysis('Params','TDR','ToggleES','OFF');
U.SetRhoDmnu2TimeEta('mnu',1,'Time',365.242*24*3600,'Eta',1);%,'RhoDMode','EffMass','RhoD',66.5e-6);
fig2=figure(2);
pT=semilogy(T.R.Te-T.R.Q,T.R.TBDDS_R,'LineWidth',2,'Color','BLue');
hold on;
pU=plot(U.R.Te-U.R.Q,U.R.TBDDS_R,'LineWidth',2,'Color','Red');
ylim([1e-15 1e-10]);
%yticks([1e-20 1e-10]);
xlabel(sprintf('{\\itE} - {\\itE}_0 (eV)'),'FontSize',12);
ylabel('Rate per 0.1 eV per second','FontSize',12);
lh1 = legend([pT pU],'Full FSD','Gaussian fit of ground state','Box','off','location','northwest');
set(lh1,'FontSize',12);
PrettyFigureFormat;
hold off;


Eta = 1e6./(T.R.NormFactorTBDDS_R*T.R.WGTS_CD_MolPerCm2/5e17./T.R.TimeSec.*365.242.*24.*3600);
T.SetRhoDmnu2TimeEta('mnu',1,'Time',365.242*24*3600,'Eta',Eta);

T.SetRhoDmnu2TimeEta('mnu',1,'Time',365.242*24*3600,'Eta',(6.4915/5)*1.7931e11/0.5717,'RhoDMode','EffMass','RhoD',66.5e-6);
T.TotalCountsStudy('Plot','OFF');

for i=1:numel(T.R.RF(1,:))
    for j=1:numel(T.R.RF(:,1))
        if T.R.RF(j,i)>0.398
            T.R.RF(j,i)=0.3983;
        end
    end
end
T.R.ComputeTBDIS;

fig2=figure(2);
hH0 = semilogy(H0(:,1)-18574,H0(:,2),'-s','LineWidth',2,'Color','Green','LineStyle','--','MarkerFaceColor','Green');
hold on;
hH1 = semilogy(H1(:,1)-18574,H1(:,2),'-s','LineWidth',2,'Color','Red','LineStyle','--','MarkerFaceColor','Red');
TBDIS_R = T.R.TBDIS./T.R.qUfrac./T.R.TimeSec;
TBDIS_B = (T.R.TBDIS-T.R.TBDIS_R)./T.R.qUfrac./T.R.TimeSec;
           
e1 = T.R.qU(T.R.qU(:,1)>(T.R.Q-310),:)-18574;
tmpis1 = TBDIS_B(T.R.qU(:,1)>(T.R.Q-310),:);
tmpis2 = TBDIS_R(T.R.qU(:,1)>(T.R.Q-310),:);

hT1 = semilogy((e1),tmpis1,'-s','LineWidth',2,'Color','Black','LineStyle','-','MarkerFaceColor','Black');
hT2 = semilogy((e1),tmpis2,'-s','LineWidth',2,'Color','Blue','LineStyle','-','MarkerFaceColor','Blue');

grid off;
%ylim([1e-20,1e2]);
xlabel('Retarding potential - 18574 (eV)','FontSize',12);
str = sprintf('Event Rate (cps)');
ylabel(str,'FontSize',14);
strT1 = sprintf('0 capture events/yr');
strT2 = sprintf('1e6 capture events/yr');
strH0 = sprintf('Formaggio 0 capture events/yr');
strH1 = sprintf('Formaggio 1e6 capture events/yr');
lh1 = legend([hT1 hT2 hH0 hH1],strT1,strT2,strH0,strH1,'Box','off');
set(lh1,'FontSize',12);
PrettyFigureFormat;
hold off;

fig3=figure(3);
p=plot(e1,(H0(:,2)-tmpis1)./tmpis1,'-s','LineWidth',2,'Color','Black','LineStyle','-');
grid on;
xlabel('E-E_0 (eV)','FontSize',12);
ylabel('Residual','FontSize',14);
str = sprintf('Residuals for 0 capture events, column density 5e17');
lh1 = legend([p], str,'Box','off');
set(lh1,'FontSize',12);
PrettyFigureFormat;


T.SetRhoDmnu2TimeEta('mnu',1,'Time',365.242*24*3600,'Eta',1.7931e11/0.5717*0.4,'RhoDMode','EffMass','RhoD',66.5e-6);
T.TotalCountsStudy('Plot','OFF');

for i=1:numel(T.R.RF(1,:))
    for j=1:numel(T.R.RF(:,1))
        if T.R.RF(j,i)>0.398
            T.R.RF(j,i)=0.3983;
        end
    end
end
T.R.ComputeTBDIS;

fig4=figure(4);
hH0 = semilogy(H0(:,1)-18574,H0(:,2),'-s','LineWidth',2,'Color','Green','LineStyle','--');
hold on;
hH1 = semilogy(H1(:,1)-18574,H1(:,2),'-s','LineWidth',2,'Color','Red','LineStyle','--');
TBDIS_R = T.R.TBDIS./T.R.qUfrac./T.R.TimeSec;
TBDIS_B = (T.R.TBDIS-T.R.TBDIS_R)./T.R.qUfrac./T.R.TimeSec;
           
e1 = T.R.qU(T.R.qU(:,1)>(T.R.Q-310),:)-18574;
tmpis1 = TBDIS_B(T.R.qU(:,1)>(T.R.Q-310),:);
tmpis2 = TBDIS_R(T.R.qU(:,1)>(T.R.Q-310),:);

hT1 = semilogy((e1),tmpis1,'-s','LineWidth',2,'Color','Black','LineStyle','-');
hT2 = semilogy((e1),tmpis2,'-s','LineWidth',2,'Color','Blue','LineStyle','-');

grid off;
%ylim([1e-20,1e2]);
xlabel('Retarding potential - 18574 (eV)','FontSize',12);
str = sprintf('Event Rate (cps)');
ylabel(str,'FontSize',14);
strT1 = sprintf('0 capture events/yr');
strT2 = sprintf('1e6 capture events/yr');
strH0 = sprintf('Formaggio 0 capture events/yr');
strH1 = sprintf('Formaggio 1e6 capture events/yr');
lh1 = legend([hT1 hT2 hH0 hH1],strT1,strT2,strH0,strH1,'Box','off');
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