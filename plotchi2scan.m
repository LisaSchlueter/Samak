load('RelicChi2Scan_TDR1.mat');

plot((0:(Netabins-1))*((etafactor*10^(etarange))/(Netabins-1)),Chi2,'LineWidth',2);
xlabel('\eta','FontSize',12);
ylabel('\chi^2','FontSize',12);
PrettyFigureFormat;