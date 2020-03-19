% plot overlay twin data for illustration

RunNr = 57089;

E0 = knm2FS_GetE0Twins('SanityPlot','OFF');
savedirT = [getenv('SamakPath'),'tritium-data/mat/TwinKnm2/'];
savedirD = [getenv('SamakPath'),'tritium-data/mat/Knm2/'];
savenameT = sprintf('%sTwin%.0f_meanE0%.2feV.mat',savedirT,RunNr,mean(E0));
savenameD = sprintf('%s%.0f.mat',savedirD,RunNr);

dT = importdata(savenameT);
dD = importdata(savenameD);

PixList = GetPixList('Knm2');
%%
f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
reNorm = sum(mean(dD.qUfrac,2)); % renorm, because of new qUfrac def.
pTwin = plot(mean(dT.qU(:,PixList),2)-18574,reNorm*sum(dT.TBDIS(:,PixList),2)./...
    mean(dT.qUfrac(:,PixList).*dT.TimeSec(PixList)',2),'LineWidth',2,'Color',rgb('Orange'));
hold on;
Rated = sum(dD.TBDIS(:,PixList),2)./mean(dD.qUfrac(:,PixList).*dD.TimeSec(PixList)',2);
RatedErr = sqrt(sum(dD.TBDIS(:,PixList),2))./mean(dD.qUfrac(:,PixList).*dD.TimeSec(PixList)',2);
pData = errorbar(mean(dD.qU(:,PixList),2)-18574,Rated,RatedErr,...
    'o','Color',rgb('DodgerBlue'),'MarkerFaceColor',rgb('DodgerBlue'),'CapSize',0,'LineWidth',2);
xlabel('Retarding energy - 18574 (eV)')
ylabel('Rate (cps)');
PrettyFigureFormat('FontSize',24)
xlim([-40 10]);
ylim([0.1,100])
leg = legend([pData,pTwin],'KATRIN Data','KATRIN MC Twin');
legend boxoff
t = title(sprintf('Run %.0f',RunNr),'FontWeight','normal');
set(gca,'YScale','log');

plotdir = [getenv('SamakPath'),'knm2ana/knm2_FigureSkating/plots/'];
plotname = sprintf('%sKnm2DataTwin%.0f.pdf',plotdir,RunNr);
export_fig(f1,plotname);
fprintf('save plot to %s \n',plotname);
