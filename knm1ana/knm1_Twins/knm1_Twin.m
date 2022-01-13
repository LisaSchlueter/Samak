
RunNr = 51883; % random run number
DataDir = [getenv('SamakPath'),'tritium-data/mat/'];
PixList = GetPixList('Knm1');

datafile = sprintf('%sKnm1/%.0f.mat',DataDir,RunNr);
data = importdata(datafile);
data_qU   = mean(data.qU(13:end,PixList),2);
data_rate = sum(data.TBDIS(13:end,PixList)')'./(mean(data.qUfrac(13:end,PixList),2).*data.TimeSec(1));
data_rate_err = sqrt(sum(data.TBDIS(13:end,PixList)')')./(mean(data.qUfrac(13:end,PixList),2).*data.TimeSec(1));

twinfile = sprintf('%sTwinKnm1/Twin%.0f_E018573.73eV.mat',DataDir,RunNr);%_BkgPtSlope-2.2muCpsPerS
twin = importdata(twinfile);
twin_qU = mean(twin.qU(13:end,PixList),2);
twin_rate = sum(twin.TBDIS(13:end,PixList)')'./(mean(twin.qUfrac(13:end,PixList),2).*twin.TimeSec(1));

%% 
close all
f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.4]);
pt = plot(twin_qU,twin_rate,'-','LineWidth',2,'Color',rgb('MediumSpringGreen'));
hold on;
ed = errorbar(data_qU,data_rate,data_rate_err,'k.','CapSize',0,'LineWidth',1.5,'MarkerSize',10);
ax = gca;
set(gca,'YScale','log');
PrettyFigureFormat('FontSize',20);
xlabel('Retarding energy (eV)');
ylabel('Rate (cps)');
ylim([0.1 40])
xlim(18574+[-42, + 50])
xticks([18535:20:18615]);
ax.XAxis.Exponent = 0;
leg = legend([ed,pt],'Data scan','MC Twin scan','Location','northeast');
legend box off
leg.FontSize = get(gca,'FontSize');
%%
pltdir = [getenv('SamakPath'),'knm1ana/knm1_Twins/plots/'];
pltname = sprintf('%sknm1_Twin_%.0f.pdf',pltdir,RunNr);
MakeDir(pltdir);
export_fig(pltname);

