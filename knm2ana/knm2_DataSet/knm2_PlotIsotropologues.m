% get Data
savedir = [getenv('SamakPath'),'knm2ana/knm2_PngBkg/results/'];
savename = sprintf('%sknm2ubfinal_Fit_Bpng-%.1fmucpsPers_%s_%.0feV_%s_%s_%s_%s_SysBudget40.mat',...
    savedir,3,'Real',40,'mNuE0BkgNorm','chi2CMShape','StackPixel','KNM2_0p1eV');
d = importdata(savename);

R = d.A;
close all

ytt = 1e2*R.SingleRunData.WGTS_MolFrac_TT;
ydt = 1e2*R.SingleRunData.WGTS_MolFrac_DT;
yht = 1e2*R.SingleRunData.WGTS_MolFrac_HT;
yt = 1e2*(R.SingleRunData.WGTS_MolFrac_TT+0.5.*R.SingleRunData.WGTS_MolFrac_HT+0.5.*R.SingleRunData.WGTS_MolFrac_DT);

LiveTime = hours(R.SingleRunData.StartTimeStamp-R.SingleRunData.StartTimeStamp(1));
x = LiveTime;
xmin= min(x)-max(x)*0.02;
xmax = max(x)*1.02;


fig88 = figure('Renderer','painters');
set(fig88,'units','normalized','pos',[1.1, 1.1,0.8,1.2]);

% put label once a week
Date = R.SingleRunData.StartTimeStamp;
TimeWeek = 7*24;
nWeek = floor(max(LiveTime)/TimeWeek);
StartIdx = 24;
myXticks = LiveTime(StartIdx)+(0:1:nWeek).*TimeWeek;
myXdates = Date(StartIdx)+(0:1:nWeek).*7;

%% tritium purity
% 
s1 =  subplot(4,6,1:5);
lt = plot(linspace(xmin,xmax,10),mean(yt)*ones(10,1),':','Color',rgb('DodgerBlue'),'LineWidth',2);
hold on;
pt = plot(x,yt,'.','MarkerSize',20,'Color',lt.Color);
PrettyFigureFormat('FontSize',LocalFontSize);
xticks(myXticks);
xticklabels('');
set(gca,'XMinorTick','off');
xlim([xmin,xmax]);
axt1 = gca;
ylabel(sprintf('%s_{T} (%%)',char(949)))
ylim(1e2*[0.978,0.996])
myYlim = ylim; 
leg = legend([pt,lt],'KNM2 scan-wise values',sprintf('\\langle%s_{T}\\rangle = %.2f %%',char(949),mean(yt)),...
    'Location','northwest');
legend box off

s2 =  subplot(4,6,6);
axt2 = gca;
h1 = histogram(yt,'FaceColor',lt.Color,'EdgeColor',lt.Color);
h1.Orientation='horizontal';
h1.BinWidth = 0.05;
ylim(myYlim);
PrettyFigureFormat('FontSize',LocalFontSize);
yticklabels('');
xticklabels('');
set(gca,'XMinorTick','off');
linkaxes([s1,s2],'y');
%% tt 
s3 =  subplot(4,6,7:11);
ltt = plot(linspace(xmin,xmax,10),mean(ytt)*ones(10,1),':','Color',rgb('FireBrick'),'LineWidth',2);
hold on;
ptt = plot(x,ytt,'.','MarkerSize',20,'Color',ltt.Color);
PrettyFigureFormat('FontSize',LocalFontSize);
xticks(myXticks);
xticklabels('');
set(gca,'XMinorTick','off');
xlim([xmin,xmax]);
axtt1 = gca;
ylabel(sprintf('c_{T_2} (%%)'));%,char(949)))
ylim(1e2*[0.961,0.992])
myYlim = ylim; 
leg = legend([ptt,ltt],'KNM2 scan-wise values',sprintf('\\langlec_{T_2}\\rangle = %.2f %%',mean(ytt)),...
    'Location','northwest');
legend box off

s4 =  subplot(4,6,12);
axtt2 = gca;
h2 = histogram(ytt,'FaceColor',ltt.Color,'EdgeColor',ltt.Color);
h2.Orientation='horizontal';
h2.BinWidth = 0.072;
ylim(myYlim);

PrettyFigureFormat('FontSize',LocalFontSize);
yticklabels('');
xticklabels('');
set(gca,'XMinorTick','off');
linkaxes([s3,s4],'y');

%% HT
s5 =  subplot(4,6,13:17);
lht = plot(linspace(xmin,xmax,10),mean(yht)*ones(10,1),':','Color',rgb('Orange'),'LineWidth',2);
hold on;
pht = plot(x,yht,'.','MarkerSize',20,'Color',lht.Color);
PrettyFigureFormat('FontSize',LocalFontSize);
xticks(myXticks);
xticklabels('');
set(gca,'XMinorTick','off');
xlim([xmin,xmax]);
axht1 = gca;
ylabel(sprintf('c_{HT} (%%)'));%,char(949)))
ylim(1e2*[0.016,0.034])
myYlim = ylim;
leg = legend([pht,lht],'KNM2 scan-wise values',sprintf('\\langlec_{HT}\\rangle = %.2f %%',mean(yht)),...
    'Location','northwest');
legend box off

s6 =  subplot(4,6,18);
axht2 = gca;
h3 = histogram(yht,'FaceColor',lht.Color,'EdgeColor',lht.Color);
h3.Orientation='horizontal';
h3.BinWidth = 0.05;
ylim(myYlim);
PrettyFigureFormat('FontSize',LocalFontSize);
yticklabels('');
xticklabels('');
set(gca,'XMinorTick','off');
linkaxes([s5,s6],'y');

%% DT
s7 =  subplot(4,6,19:23);
ldt = plot(linspace(xmin,xmax,10),mean(ydt)*ones(10,1),':','Color',rgb('ForestGreen'),'LineWidth',2);
hold on;
pdt = plot(x,ydt,'.','MarkerSize',20,'Color',ldt.Color);
PrettyFigureFormat('FontSize',LocalFontSize);
xlabel(sprintf('Date in %s',datestr(R.SingleRunData.StartTimeStamp(1),'yyyy')));
xticks(myXticks);
xticklabels(datestr(myXdates,'mmm, dd'));
set(gca,'XMinorTick','off');
xlim([xmin,xmax]);
axdt1 = gca;
ylabel(sprintf('c_{DT} (%%)'));%,char(949)))
ylim(1e2*[0.0021,0.0044])
myYlim = ylim;
leg = legend([pdt,ldt],'KNM2 scan-wise values',sprintf('\\langlec_{DT}\\rangle = %.2f %%',mean(ydt)),...
    'Location','northwest');
legend box off

s8 =  subplot(4,6,24);
axdt2 = gca;
h4 = histogram(ydt,'FaceColor',ldt.Color,'EdgeColor',ldt.Color);
h4.Orientation='horizontal';
h4.BinWidth = 0.007;
ylim(myYlim);
PrettyFigureFormat('FontSize',LocalFontSize);
yticklabels('');
xlabel('Entries');
set(gca,'XMinorTick','off');
linkaxes([s5,s6],'y');

linkaxes([s1,s3,s5,s7],'x');
linkaxes([s2,s4,s6,s8],'x');
xlim([0 65])
%% change positions
outerpos1 = axtt1.OuterPosition;
ti = axtt1.TightInset;
left = outerpos1(1) + ti(1);
bottom1 = outerpos1(2) + ti(2);
ax_width = outerpos1(3) - ti(1);% - ti(3);
ax_height = outerpos1(4)+0.03;% - ti(2) - ti(4);
outerpos2 = axht1.OuterPosition;
bottom2 = outerpos2(2) + ti(2);
outerpos3 = axdt1.OuterPosition;
bottom3 = outerpos3(2) + ti(2);
%%
axt1.Position = [left-0.03 bottom1+0.255 ax_width+0.02 ax_height];
axtt1.Position = [left-0.03 bottom1+0.05 ax_width+0.02 ax_height];
axht1.Position = [left-0.03 bottom2+0.062 ax_width+0.02 ax_height];
axdt1.Position = [left-0.03 bottom3+0.12 ax_width+0.02 ax_height];

%%
%adjust axes to be adjacent to first subplot
axt2.Position(1) = axt1.Position(1)+axt1.Position(3).*1.01;%left
axt2.Position(2) = axt1.Position(2); % bottom
axt2.Position(4) = axt1.Position(4); % height
axt2.Position(3) = 0.13;
%adjust axes to be adjacent to first subplot
axtt2.Position(1) = axtt1.Position(1)+axtt1.Position(3).*1.01;%left
axtt2.Position(2) = axtt1.Position(2); % bottom
axtt2.Position(4) = axtt1.Position(4); % height
axtt2.Position(3) = 0.13;
%adjust axes to be adjacent to first subplot
axht2.Position(1) = axht1.Position(1)+axht1.Position(3).*1.01;%left
axht2.Position(2) = axht1.Position(2); % bottom
axht2.Position(4) = axht1.Position(4); % height
axht2.Position(3) = 0.13;
%adjust axes to be adjacent to first subplot
axdt2.Position(1) = axdt1.Position(1)+axdt1.Position(3).*1.01;%left
axdt2.Position(2) = axdt1.Position(2); % bottom
axdt2.Position(4) = axdt1.Position(4); % height
axdt2.Position(3) = 0.13;

%%
pltdir = [getenv('SamakPath'),'knm2ana/knm2_DataSet/plots/'];
pltfile = sprintf('%sknm2_Isotopologues.pdf',pltdir);
export_fig(pltfile);