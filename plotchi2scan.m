load('RelicChi2Scan_TDR4.mat');

%Chi2(3)=0.16;
fig1=figure(1);
plot((0:(Netabins-1))*((etafactor*10^(etarange))/(Netabins-1)),Chi2,'LineWidth',2);
xlabel('\eta','FontSize',12);
ylabel('\chi^2','FontSize',12);
PrettyFigureFormat;

% t=ref_RelicNuBkg_DesignReport('eta_i',2.81e9);
% t.ComputeNormFactorTBDDS;
% t.ComputeTBDDS;
% t.ComputeTBDIS;
% 
% TBDIS_R = t.TBDIS./t.qUfrac./t.TimeSec;
% TBDIS_B = (t.TBDIS-t.TBDIS_R)./t.qUfrac./t.TimeSec;
%            
% e1 = t.qU(t.qU(:,1)>(t.Q-310),:)-18575;
% tmpisB = TBDIS_B(t.qU(:,1)>(t.Q-310),:);
% tmpisR = TBDIS_R(t.qU(:,1)>(t.Q-310),:);
% fig2=figure(2);
% p1=plot(e1,tmpisR./tmpisB,'-s','LineWidth',2,'LineStyle','-');
% hold on;
% 
% t=ref_RelicNuBkg_DesignReport('eta_i',2.81e9,'BKG_RateAllFPDSec',0.05);
% t.ComputeNormFactorTBDDS;
% t.ComputeTBDDS;
% t.ComputeTBDIS;
% 
% TBDIS_R = t.TBDIS./t.qUfrac./t.TimeSec;
% TBDIS_B = (t.TBDIS-t.TBDIS_R)./t.qUfrac./t.TimeSec;
%            
% e1 = t.qU(t.qU(:,1)>(t.Q-310),:)-18575;
% tmpisB = TBDIS_B(t.qU(:,1)>(t.Q-310),:);
% tmpisR = TBDIS_R(t.qU(:,1)>(t.Q-310),:);
% p2=plot(e1,tmpisR./tmpisB,'-s','LineWidth',2,'LineStyle','-');
% 
% t=ref_RelicNuBkg_DesignReport('eta_i',2.81e9,'BKG_RateAllFPDSec',0.1);
% t.ComputeNormFactorTBDDS;
% t.ComputeTBDDS;
% t.ComputeTBDIS;
% 
% TBDIS_R = t.TBDIS./t.qUfrac./t.TimeSec;
% TBDIS_B = (t.TBDIS-t.TBDIS_R)./t.qUfrac./t.TimeSec;
%            
% e1 = t.qU(t.qU(:,1)>(t.Q-310),:)-18575;
% tmpisB = TBDIS_B(t.qU(:,1)>(t.Q-310),:);
% tmpisR = TBDIS_R(t.qU(:,1)>(t.Q-310),:);
% p3=plot(e1,tmpisR./tmpisB,'-s','LineWidth',2,'LineStyle','-');
% 
% t=ref_RelicNuBkg_DesignReport('eta_i',2.81e9,'BKG_RateAllFPDSec',0.15);
% t.ComputeNormFactorTBDDS;
% t.ComputeTBDDS;
% t.ComputeTBDIS;
% 
% TBDIS_R = t.TBDIS./t.qUfrac./t.TimeSec;
% TBDIS_B = (t.TBDIS-t.TBDIS_R)./t.qUfrac./t.TimeSec;
%            
% e1 = t.qU(t.qU(:,1)>(t.Q-310),:)-18575;
% tmpisB = TBDIS_B(t.qU(:,1)>(t.Q-310),:);
% tmpisR = TBDIS_R(t.qU(:,1)>(t.Q-310),:);
% p4=plot(e1,tmpisR./tmpisB,'-s','LineWidth',2,'LineStyle','-');
% 
% t=ref_RelicNuBkg_DesignReport('eta_i',2.81e9,'BKG_RateAllFPDSec',0.2);
% t.ComputeNormFactorTBDDS;
% t.ComputeTBDDS;
% t.ComputeTBDIS;
% 
% TBDIS_R = t.TBDIS./t.qUfrac./t.TimeSec;
% TBDIS_B = (t.TBDIS-t.TBDIS_R)./t.qUfrac./t.TimeSec;
%            
% e1 = t.qU(t.qU(:,1)>(t.Q-310),:)-18575;
% tmpisB = TBDIS_B(t.qU(:,1)>(t.Q-310),:);
% tmpisR = TBDIS_R(t.qU(:,1)>(t.Q-310),:);
% p5=plot(e1,tmpisR./tmpisB,'-s','LineWidth',2,'LineStyle','-');
% xlabel('E-E_{0} (eV)','FontSize',12);
% ylabel('S/B');
% %ylim([1,1.008]);
% lh1=legend([p1 p2 p3 p4 p5],'10 mcps','50 mcps','100 mcps','150 mcps','200 mcps');
% set(lh1,'FontSize',12);
% PrettyFigureFormat;
% hold off;