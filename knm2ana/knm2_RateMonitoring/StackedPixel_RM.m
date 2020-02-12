%
% KNM2 - -300V Rate Monitor FPD
%  
% Stakced-pixel Evolution 
%
% Last Modified: 17/10/2019
% T. Lasserre
% 
%

%% Read Data
DataType = 'Real';
RunList = 'KNM2_RW3';
FSDFlag = 'BlindingKNM2';
ELossFlag = 'KatrinT2';
AnaFlag = 'StackPixel'; % uniform FPD
RunAnaArg = {'RunList',RunList,'DataType',DataType,...
            'FSDFlag',FSDFlag,'ELossFlag',ELossFlag,'AnaFlag',AnaFlag};
A = MultiRunAnalysis(RunAnaArg{:});

%% Slow Control Data
p1 =(A.SingleRunData.WGTS_MolFrac_TT'+0.5*A.SingleRunData.WGTS_MolFrac_HT'+0.5*A.SingleRunData.WGTS_MolFrac_DT')./mean((A.SingleRunData.WGTS_MolFrac_TT'+0.5*A.SingleRunData.WGTS_MolFrac_HT'+0.5*A.SingleRunData.WGTS_MolFrac_DT')).*A.SingleRunData.WGTS_CD_MolPerCm2'./mean(A.SingleRunData.WGTS_CD_MolPerCm2');
p2 = mean(A.SingleRunData.qU_RM,1); p2=p2-mean(p2);           

%% Stacked Pixel Data
count  = sum(A.SingleRunData.TBDIS_RM(A.PixList,:),1);
sstime = mean(A.SingleRunData.qUfrac_RM,1).*A.SingleRunData.TimeSec;
rate   = count./sstime;
cf     = A.RMRateErosCorrectionqUActivity;
%cf     = p1';

%% Stacked Pixel Histogram - Raw Rate
figure(1)
subplot(2,1,1)
histogram(rate);
PrettyFigureFormat
xlabel('rate (cps)');
ylabel('beta-scans');
title(sprintf('Raw rate at qU = -%.1f V',mean(mean(A.SingleRunData.qU_RM))));
subplot(2,1,2)
histogram(rate.*cf);
PrettyFigureFormat
xlabel('rate (cps)');
ylabel('beta-scans');
title(sprintf('Corrected rate at qU = -%.1f V',mean(mean(A.SingleRunData.qU_RM))));


%% Stacked Pixel Evolution - Raw Rate - Corrected Rate
figure(11)
subplot(2,1,1)
plot(p1,'s-','Color',rgb('DarkBlue'),'LineWidth',1,'MarkerSize',4,'markerfacecolor','Black');
ylabel('relative activity');
xlabel('beta-scan');
PrettyFigureFormat
subplot(2,1,2)
plot(p2,'s--','Color',rgb('DarkGreen'),'LineWidth',1,'MarkerSize',8,'markerfacecolor',rgb('DarkGreen'));
ylabel('relative qU');
xlabel('beta-scan');
PrettyFigureFormat
title(sprintf('Slow Control Parameters'));


%% Stacked Pixel Evolution : Raw Rate - Corrected Rate
figure(2)
subplot(2,1,1)
plot((cf-1).*rate,'s-','Color',rgb('DarkBlue'),'LineWidth',1,'MarkerSize',4,'markerfacecolor','Black');
ylabel('rate variation (cps)');
xlabel('beta-scan');
PrettyFigureFormat
subplot(2,1,2)
plot(rate,'s--','Color',rgb('DarkGreen'),'LineWidth',1,'MarkerSize',8,'markerfacecolor',rgb('DarkGreen'));
hold on
plot(rate.*cf,'s--','Color',rgb('IndianRed'),'LineWidth',1,'MarkerSize',8,'markerfacecolor',rgb('IndianRed'));
hold off
ylabel('rate (cps)');
xlabel('beta-scan');
PrettyFigureFormat
title(sprintf('Raw rate at qU = -%.1f V',mean(mean(A.SingleRunData.qU_RM))));

%% Stacked Pixel: conversion to counts
count_norm     = rate .* mean(sstime);
corrcount_norm = count_norm .* cf;
figure(3)
subplot(2,1,1)
histogram(count_norm);
PrettyFigureFormat
xlabel('counts');
ylabel('beta-scans');
title(sprintf('Raw rate at qU = -%.1f V',mean(mean(A.SingleRunData.qU_RM))));
subplot(2,1,2)
histogram(corrcount_norm.*cf);
PrettyFigureFormat
xlabel('corrected counts');
ylabel('beta-scans');
title(sprintf('Corrected rate at qU = -%.1f V',mean(mean(A.SingleRunData.qU_RM))));

%% Stacked Pixel: Non-Poissonian Component?
pdG = fitdist(corrcount_norm','Normal');
pdN = fitdist(corrcount_norm','poisson');

myMainTitle = sprintf('KATRIN - %s - FPD Non-Poisson Rate 300eV below Endpoint',RunList);
maintitle   = myMainTitle;
savefile    = sprintf('plots/test.png');
fig1      = figure('Name',sprintf('KATRIN - %s Scanwise Background',RunList),...
    'NumberTitle','off','rend','painters','pos',[10 10 1400 1000]);
a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
a.FontSize=24;a.FontWeight='bold';
h = histogram(corrcount_norm,15,'Normalization','pdf',...
    'FaceColor',rgb('DodgerBlue'),'LineWidth',2,'FaceAlpha',0.7);
xlabel(sprintf('counts in %.2f sec',mean(sstime)));
ylabel('Frequency');
PrettyFigureFormat
% Gaussian PDF
pdfG = @(b) 1/sqrt(2*pi*pdG.sigma.^2) .* exp(-(b-pdG.mu).^2/(2*pdG.sigma.^2));
% Poisson PDF
pdfP = @(b) 1/sqrt(2*pi*pdN.lambda) .* exp(-(b-pdN.lambda).^2/(2*pdN.lambda));
hold on
%b = (h.BinEdges(1:end-1)+h.BinWidth/2);
b = linspace(min(corrcount_norm),max(corrcount_norm),100);
g=plot(b,pdfG(b),'Color',rgb('IndianRed'),'LineWidth',4);
p=plot(b,pdfP(b),'Color',rgb('GoldenRod'),'LineWidth',4);
leg=legend([h g p],...
    sprintf('subscans'),...
    sprintf('Gaussian \\sigma=%.2f counts',...
    (pdG.sigma)),sprintf('Poisson \\sigma=%.2f counts',sqrt(pdN.lambda)),...
    'location','northwest');
leg.Color = 'none'; legend boxoff;
hold off
PrettyFigureFormat
set(gca,'FontSize',24);
export_fig(gcf,savefile,'-m3');
disp(pdG.sigma/sqrt(pdN.lambda));

%% Rate Evolution
myMainTitle = sprintf('KATRIN - %s - FPD Rate Evolution 300eV below Endpoint',RunList);
maintitle   = myMainTitle;
savefile    = sprintf('plots/test.png');
fig2      = figure('Name',sprintf('KATRIN - %s Scanwise Background',RunList),...
    'NumberTitle','off','rend','painters','pos',[10 10 1400 1000]);
a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
a.FontSize=24;a.FontWeight='bold';

subplot(3,1,1)
plot(A.SingleRunData.StartTimeStamp,p1,'s-','Color',rgb('DarkBlue'),'LineWidth',1,'MarkerSize',12,'markerfacecolor','Black');
ylabel('relative activity');
xlabel('Scan Start Time');
PrettyFigureFormat

subplot(3,1,2)
plot(A.SingleRunData.StartTimeStamp,p2,'s--','Color',rgb('DarkGreen'),'LineWidth',1,'MarkerSize',12,'markerfacecolor',rgb('DarkGreen'));
ylabel('relative qU');
xlabel('Scan Start Time');
PrettyFigureFormat

subplot(3,1,3)
hnc=plot(A.SingleRunData.StartTimeStamp,rate,'s--','Color',rgb('DarkGreen'),'LineWidth',1,'MarkerSize',8,'markerfacecolor',rgb('DarkGreen'));
hold on
hc=plot(A.SingleRunData.StartTimeStamp,rate.*cf,'s--','Color',rgb('IndianRed'),'LineWidth',1,'MarkerSize',12,'markerfacecolor',rgb('IndianRed'));
hold off
ylabel('rate (cps)');
xlabel('Scan Start Time');
leg=legend([hnc hc],...
    sprintf('before correction'),...
    sprintf('after correction'),...
    'location','northwest');
leg.Color = 'none'; legend boxoff;
PrettyFigureFormat

