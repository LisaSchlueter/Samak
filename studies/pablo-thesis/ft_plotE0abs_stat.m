addpath(genpath('../../../Samak2.0'));
close all
%BuildTable(1);

% load('BuildTableWorkspace2')

%(TotUncVec,StaUncVec,E_totalSys);
yn = 1:12;
% E_totalStat = flipud(E_totalStat);
% StaUncVec = flipud(StaUncVec);

yana = ["Unseg. Av. Runs Short","Unseg. Av. Runs Med.",...
    "Unseg. Av. Runs Long","Unseg. Stack. Runs Short",...
    "Unseg. Stack. Runs Med.","Unseg. Stack. Runs Long",...
    "Ring Short","Ring Med.","Single Pix. Short","Single Pix. Med.",...
    "Multi Pix. Short","Multi Pix. Med."]';

e0x = string(18573.3:0.4:18574.7);

yana = flipud(yana);
fige0 = figure(659);

pos = get(fige0,'position');
set(fige0,'position',[pos(1:2)/3.5 pos(3:4)*1.5])
axes('YAxisLocation','right')
workFun = 1;

hold on
h1 = plot(E_totalStat,yn,'sk','MarkerFaceColor','black');
errorb(E_totalStat,yn,StaUncVec,'horizontal','color','black','barwidth',0.125);
htheo = plot(18573.7*ones(1,length(E_totalStat)+7),0:18,'Color','blue','LineStyle','--','LineWidth',1);
hliml = plot((18573.7-workFun/2)*ones(1,length(E_totalStat)+7),0:18,'Color',rgb('CadetBlue'),'LineStyle','-','LineWidth',0.5);
hlimu = plot((18573.7 + workFun/2)*ones(1,length(E_totalStat)+7),0:18,'Color',rgb('CadetBlue'),'LineStyle','-','LineWidth',0.5);
rectangle('Position',[18573.7-workFun/2,0,workFun,15.8],'Facecolor',[rgb('CadetBlue'), 0.5],...
    'EdgeColor','none');
hrec=fill([1 1 0 0],[0 1 1 0],rgb('CadetBlue'));
hold off

axis([18572.7,18574.7,0,15.8])
legend([h1,htheo,hrec],{'Statistical Uncertainty','Expected Value',...
    'Work Function Uncertainty'},'location','northeast')
legend('boxoff')
PrettyFigureFormat;
yticks(yn)
yticklabels(yana);

%ytickangle(90);
xticks(18572.9:0.4:18574.7)
xticklabels(["(eV)",e0x]);

title('Effective endpoint using only statistics in all analysis types')

export_fig('plots/E0SummaryAbsoluteStat.pdf','-pdf')