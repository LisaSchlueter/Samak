AnaFlag   = 'StackPixel';

savedir = [getenv('SamakPath'),'knm2ana/knm2_RateMonitoring/results/'];
savename = sprintf('%sknm2_RateMonitoring_GetCounts_%s.mat',savedir,AnaFlag);
load(savename);

%%
x = LiveTime;
xmin= min(x)-max(x)*0.02;
xmax = max(x)*1.02;

% put label once a week
TimeWeek = 7*24;
nWeek = floor(max(LiveTime)/TimeWeek);
StartIdx = 24;
myXticks = LiveTime(StartIdx)+(0:1:nWeek).*TimeWeek;
myXdates = StartTimeStamp(StartIdx)+(0:1:nWeek).*7;

RatesErr = sqrt(Rates.*TimeSec_ScanStep)./TimeSec_ScanStep;
RatesCorrErr = sqrt(RatesCorr.*TimeSec_ScanStep)./TimeSec_ScanStep;

fig88 = figure('Renderer','painters');
set(fig88,'units','normalized','pos',[1.1, 1.1,0.8,0.6]);

RateArg = {'.','MarkerSize',12,'LineWidth',1.5,'CapSize',0};
s1 = subplot(5,1,[1:3]);
% uncorrected
errorbar(x(Idx_Rw1),1e-04.*Rates(Idx_Rw1),1e-04.*RatesErr(Idx_Rw1),RateArg{:},'Color',rgb('DimGray'));
hold on;
errorbar(x(Idx_Rw2),1e-04.*Rates(Idx_Rw2),1e-04.*RatesErr(Idx_Rw2),RateArg{:},'Color',rgb('DimGray'));
errorbar(x(Idx_Rw3),1e-04.*Rates(Idx_Rw3),1e-04.*RatesErr(Idx_Rw3),RateArg{:},'Color',rgb('DimGray'));

% corrected
pc1 = errorbar(x(Idx_Rw1),1e-04.*RatesCorr(Idx_Rw1),1e-04.*RatesCorrErr(Idx_Rw1),RateArg{:},'Color',rgb('Crimson'));
pc2 = errorbar(x(Idx_Rw2),1e-04.*RatesCorr(Idx_Rw2),1e-04.*RatesCorrErr(Idx_Rw2),RateArg{:},'Color',rgb('DodgerBlue'));
pc3 = errorbar(x(Idx_Rw3),1e-04.*RatesCorr(Idx_Rw3),1e-04.*RatesCorrErr(Idx_Rw3),RateArg{:},'Color',rgb('Orange'));

% means
MeanRate1 = mean(RatesCorr(Idx_Rw1));
MeanRate2 = mean(RatesCorr(Idx_Rw2));
MeanRate3 = mean(RatesCorr(Idx_Rw3));
pm1 = plot(x(Idx_Rw1),1e-04.*MeanRate1.*ones(sum(Idx_Rw1),1),'-','Color',pc1.Color,'LineWidth',2);
pm2 = plot(x(Idx_Rw2),1e-04.*MeanRate2.*ones(sum(Idx_Rw2),1),'-','Color',pc2.Color,'LineWidth',2);
pm3 = plot(x(Idx_Rw3),1e-04.*MeanRate3.*ones(sum(Idx_Rw3),1),'-','Color',pc3.Color,'LineWidth',2);

PrettyFigureFormat;
xticks(myXticks);
xticklabels('');
ylabel(sprintf('Rate (10^4 \\times cps)'));
set(gca,'XMinorTick','off');
xlim([xmin,xmax])
ylim([6.635 6.717]);
ax1 = gca;

leg = legend([pm1,pm2,pm3],...
    sprintf('Period 1: \\mu = %.3f \\times 10^4 cps',MeanRate1.*1e-04),...
    sprintf('Period 2: \\mu = %.3f \\times 10^4 cps',MeanRate2.*1e-04),...
    sprintf('Period 3: \\mu = %.3f \\times 10^4 cps',MeanRate3.*1e-04),...
    'Location','south');
leg.NumColumns=3;
PrettyLegendFormat(leg);
leg.FontSize = get(gca,'FontSize');
%
s2 = subplot(5,1,4);
CorrArg = {'-','LineWidth',2};
plot(linspace(min(x)-10,max(x)+10,10),zeros(10,1),':','Color',rgb('Silver'),'LineWidth',2);
hold on;
pa1 = plot(x(Idx_Rw1),1e2.*(Activity(Idx_Rw1)-mean(Activity))./mean(Activity),CorrArg{:},'Color',pc1.Color);
pa2 = plot(x(Idx_Rw2),1e2.*(Activity(Idx_Rw2)-mean(Activity))./mean(Activity),CorrArg{:},'Color',pc2.Color);
pa3 = plot(x(Idx_Rw3),1e2.*(Activity(Idx_Rw3)-mean(Activity))./mean(Activity),CorrArg{:},'Color',pc3.Color);
ylabel(sprintf('\\Delta{\\itA} (%%)'));
PrettyFigureFormat;
xticks(myXticks);
xticklabels('');
set(gca,'XMinorTick','off');
ax2 = gca;

s3 = subplot(5,1,5);
CorrArg2 = {'-.','LineWidth',2};
plot(linspace(min(x)-10,max(x)+10,10),zeros(10,1),':','Color',rgb('Silver'),'LineWidth',2);
hold on;
pa1 = plot(x(Idx_Rw1),1e3.*qUDiff(Idx_Rw1),CorrArg2{:},'Color',pc1.Color);
pa2 = plot(x(Idx_Rw2),1e3.*qUDiff(Idx_Rw2),CorrArg2{:},'Color',pc2.Color);
pa3 = plot(x(Idx_Rw3),1e3.*qUDiff(Idx_Rw3),CorrArg2{:},'Color',pc3.Color);
xlabel(sprintf('Date in %s',datestr(StartTimeStamp(1),'yyyy')));
ylabel(sprintf('\\Delta{\\itqU} (meV)'));
PrettyFigureFormat;
xticks(myXticks);
xticklabels(datestr(myXdates,'mmm, dd'));
set(gca,'XMinorTick','off');
xlim([xmin,xmax])
ax3 = gca;
ylim([-17 14]);
%
linkaxes([s1,s2,s3],'x');
xlim([xmin,xmax])

ax1.Position(4) = 0.45;
ax1.Position(2) = 0.53;
ax2.Position(2) = 0.31;
ax3.Position(2) = 0.09;

ax2.Position(4) = 0.2;
ax3.Position(4) = ax2.Position(4);

ax1.YLabel.Position(1) = -90;
ax2.YLabel.Position(1) = ax1.YLabel.Position(1)+15;
ax3.YLabel.Position(1) = ax2.YLabel.Position(1);


pltdir = [getenv('SamakPath'),'knm2ana/knm2_RateMonitoring/plots/'];
MakeDir(pltdir);
pltname = sprintf('%sknm2_RateMonitoring_GetCounts_%s.pdf',pltdir,AnaFlag);
export_fig(pltname);
fprintf('save plot to %s \n',pltname);