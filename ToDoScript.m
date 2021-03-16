PlotBestFit('DataType','Twin','Nfit',500,'Plot','OFF','Syst','ON');%,'pullFlag',[3 22 23]);

%A=RelicNuAnalysis('Params','KNM1');
%A.EtaFit('NmNuBins',5,'MaxmNu',1,'DataType','Real','Mode','SinglemnuSq','fitPar','E0 Norm Bkg eta','DeltaChi2',1,'Recompute','OFF');
%A.EtaFit('NmNuBins',5,'MaxmNu',1,'DataType','Real','Mode','SinglemnuSq','fitPar','E0 Norm Bkg eta','DeltaChi2',4,'Recompute','OFF');
%A.EtaFit('NmNuBins',5,'MaxmNu',1,'DataType','Real','Mode','SinglemnuSq','fitPar','E0 Norm Bkg eta','DeltaChi2',9,'Recompute','OFF');
%A.Chi2Scan_2D('Syst','ON','DataType','Real','DeltaEta',5e11,'Deltamnu',3,'Recompute','OFF');

%A.Chi2Twin('Recompute','OFF','CheckErrors','ON','Plot','OFF','Syst','ON','pullFlag',[1 3]);
%A.Chi2Twin('Recompute','OFF','Plot','OFF','Syst','ON','pullFlag',[1 3]);
%A.Chi2Twin('Recompute','OFF','Plot','OFF','Syst','OFF','fitPar','E0 Norm Bkg');
%fig1=figure(1);
%A.SystBreakdown('TwinBias_mnuSq',0);
%fig2=figure(2);
%A.SystBreakdown('TwinBias_mnuSq',1);
%fig3=figure(3);
%A.SystBreakdown('TwinBias_mnuSq',5);
%A.Chi2Scan_2D('Netabins',10,'Nmnubins',10,'TwinBias_mnuSq',0)
% load('./RelicNuBkg/Chi2Scans/RelicChi2Scan_Twin_BiasmnuSq1_SystON_range40_KNM1_[0 5e+11]_mNu E0 Norm Bkg_mnuSqPull.mat');
% Bkg_pull       = Bkg;
% Bkg_err_pull   = Bkg_err;
% Chi2_pull      = Chi2;
% E0_pull        = E0;
% E0_err_pull    = E0_err;
% mnuSq_pull     = mnuSq;
% mnuSq_err_pull = mnuSq_err;
% Norm_pull      = Norm;
% Norm_err_pull  = Norm_err;
% load('./RelicNuBkg/Chi2Scans/RelicChi2Scan_Twin_BiasmnuSq1_SystON_range40_KNM1_[0 5e+11]_mNu E0 Norm Bkg.mat');
% FitStyleArg = {'o','Color','k','LineWidth',1.0,'MarkerFaceColor',rgb('Black'),'MarkerSize',4,'MarkerEdgeColor',rgb('Black')};
% FitStyleArgP = {'o','Color','r','LineWidth',1.0,'MarkerFaceColor',rgb('Red'),'MarkerSize',4,'MarkerEdgeColor',rgb('Red')};
% eta=(0:(Netabins-1))*((etafactor*10^(etarange))/(Netabins-1));
% 
% fig1=figure(1);
% n=errorbar(eta,mnuSq,mnuSq_err,FitStyleArg{:},'CapSize',0);
% hold on;
% p=errorbar(eta,mnuSq_pull,mnuSq_err_pull,FitStyleArgP{:},'CapSize',0);
% fitobject  = fit(permute(eta,[2 1]),permute(mnuSq,[2 1]),'poly1');
% fitobjectP = fit(permute(eta,[2 1]),permute(mnuSq_pull,[2 1]),'poly1');
% plot(fitobject,eta,mnuSq);
% plot(fitobjectP,eta,mnuSq_pull);
% lgd = legend([n p],'No pull','Pull','Location','northwest');
% legend boxoff;
% xlabel('\eta','FontSize',12);
% ylabel('m_{\nu}^{2} (eV^{2})','FontSize',12);
% PrettyFigureFormat;
% hold off;
% 
% fig2=figure(2);
% n=errorbar(eta,E0,E0_err,FitStyleArg{:},'CapSize',0);
% hold on;
% p=errorbar(eta,E0_pull,E0_err_pull,FitStyleArgP{:},'CapSize',0);
% fitobject  = fit(permute(eta,[2 1]),permute(E0,[2 1]),'poly1');
% fitobjectP = fit(permute(eta,[2 1]),permute(E0_pull,[2 1]),'poly1');
% plot(fitobject,eta,E0);
% plot(fitobjectP,eta,E0_pull);
% lgd = legend([n p],'No pull','Pull','Location','northwest');
% legend boxoff;
% xlabel('\eta','FontSize',12);
% ylabel('E_{0}^{fit} (eV)','FontSize',12);
% PrettyFigureFormat;
% hold off;
% 
% fig3=figure(3);
% n=errorbar(eta,Norm,Norm_err,FitStyleArg{:},'CapSize',0);
% hold on;
% p=errorbar(eta,Norm_pull,Norm_err_pull,FitStyleArgP{:},'CapSize',0);
% fitobject  = fit(permute(eta,[2 1]),permute(Norm,[2 1]),'poly1');
% fitobjectP = fit(permute(eta,[2 1]),permute(Norm_pull,[2 1]),'poly1');
% plot(fitobject,eta,Norm);
% plot(fitobjectP,eta,Norm_pull);
% lgd = legend([n p],'No pull','Pull','Location','northwest');
% legend boxoff;
% xlabel('\eta','FontSize',12);
% ylabel('N','FontSize',12);
% PrettyFigureFormat;
% hold off;
% 
% fig4=figure(4);
% n=errorbar(eta,Bkg,Bkg_err,FitStyleArg{:},'CapSize',0);
% hold on;
% p=errorbar(eta,Bkg_pull,Bkg_err_pull,FitStyleArgP{:},'CapSize',0);
% fitobject  = fit(permute(eta,[2 1]),permute(Bkg,[2 1]),'poly1');
% fitobjectP = fit(permute(eta,[2 1]),permute(Bkg_pull,[2 1]),'poly1');
% plot(fitobject,eta,Bkg);
% plot(fitobjectP,eta,Bkg_pull);
% lgd = legend([n p],'No pull','Pull','Location','northwest');
% legend boxoff;
% xlabel('\eta','FontSize',12);
% ylabel('Background (cps)','FontSize',12);
% PrettyFigureFormat;
% hold off;