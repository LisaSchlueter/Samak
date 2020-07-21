load('RelicChi2Scan.mat');

plot(0.1.*1e6.*(0:Netabins-1),Chi2,'LineWidth',2);
xlabel('eta','FontSize',12);
ylabel('chi^2','FontSize',12);
PrettyFigureFormat;