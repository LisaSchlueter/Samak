% plot overlay: twin and real data for 1 run
R = RunAnalysis('RunNr',51936,'DataType','Real','exclDataStart',14);
T =RunAnalysis('RunNr',51936,'DataType','Twin','exclDataStart',14);
%%
fig6 = figure(6);
set(fig6, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.5]);
TimeSubrun = (R.RunData.qUfrac(R.exclDataStart:end).*R.RunData.TimeSec);
twin = plot(T.RunData.qU(R.exclDataStart:end)-18573.7,T.RunData.TBDIS(R.exclDataStart:end)./TimeSubrun,'-');
twin.LineWidth = 2; twin.Color = T.PlotColor;
hold on;
data = errorbar(R.RunData.qU(R.exclDataStart:end)-18573.7,R.RunData.TBDIS(R.exclDataStart:end)./TimeSubrun,...
    R.RunData.TBDISE(R.exclDataStart:end)./TimeSubrun,'--o');
data.LineWidth = 3; data.Color = rgb('RoyalBlue'); data.MarkerFaceColor = rgb('CadetBlue'); data.MarkerSize = 9;
hold off;
leg = legend([data,twin],sprintf('data run %.0f',R.RunNr),sprintf('twin run %.0f',R.RunNr)); 
legend boxoff;
PrettyFigureFormat;
xlabel('retarding potential - 18573.7 (V)');
ylabel('rate (cps)')
set(gca,'YScale','log');
ylim([min(R.RunData.TBDIS(R.exclDataStart:end)./TimeSubrun)*0.8,max(R.RunData.TBDIS(R.exclDataStart:end)./TimeSubrun)]);
set(gca,'FontSize',24);
savedir = [getenv('SamakPath'),'knm1ana/knm1Twins/plots/'];
MakeDir(savedir);
savename = sprintf('Twin%.0f_Overlay.png',R.RunNr);
print([savedir,savename],'-dpng','-r450');