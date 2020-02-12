addpath(genpath('../../../Samak2.0'));
close all
%BuildTable(1);

load('BuildTableWorkspace')

%(TotUncVec,StaUncVec,E_totalSys);
yn = 1:12;
E_totalSys = flipud(E_totalSys);
StaUncVec = flipud(StaUncVec);
TotUncVec = flipud(TotUncVec);

yana = ["UAR Short","UAR Med. (*)",...
    "UAR Long (*)","USR Short (*)",...
    "USR Med. (*)","USR Long (*)",...
    "Ring Short","Ring Med. (*)","Single Pix. Short","Single Pix. Med.",...
    "Multi Pix. Short (*)","Multi Pix. Med. (*)"]';

% yana = ["Unseg. Av. Runs Short","Unseg. Av. Runs Med. (*)",...
%     "Unseg. Av. Runs Long (*)","Unseg. Stack. Runs Short",...
%     "Unseg. Stack. Runs Med.","Unseg. Stack. Runs Long",...
%     "Ring Short (*)","Ring Med. (*)","Single Pix. Short","Single Pix. Med.",...
%     "Multi Pix. Short","Multi Pix. Med."]';

e0x = string(18573:18576);

yana = flipud(yana);
fige0 = figure(659);

pos = get(fige0,'position');
set(fige0,'position',[pos(1:2)/4 pos(3:4)*1.75])
axes('YAxisLocation','right')
workFun = 1;

hold on
%badp = [11,10,6,5]; goodp = [12,9,8,7,4,3,2,1];
badp = [11,10,6,5]; goodp = [12,9,8,7,4,3,2,1];
h1 = plot(E_totalSys,yn,'sk','MarkerFaceColor','red');
h2 = plot(E_totalSys,yn,'sk','MarkerFaceColor','black');
errorb(E_totalSys(goodp),yn(goodp),StaUncVec(goodp),'horizontal','color','red','barwidth',0.125);
errorb(E_totalSys(goodp),yn(goodp),TotUncVec(goodp),'horizontal','barwidth',0.125);
errorb(E_totalSys(badp),yn(badp),StaUncVec(badp),'horizontal','color',[255,99,71]/256,'barwidth',1);
errorb(E_totalSys(badp),yn(badp),TotUncVec(badp),'horizontal','color',[100,100,100]/256,'barwidth',1);
htheo = plot(18573.7*ones(1,length(E_totalSys)+7),0:18,'Color','blue','LineStyle','--','LineWidth',1);
% hliml = plot((E_totalSys(end)-TotUncVec(end))*ones(1,length(E_totalSys)+7),0:18,'Color',[0.5 0.5 0.5],'LineStyle','-','LineWidth',3);
% hlimu = plot((E_totalSys(end)+TotUncVec(end))*ones(1,length(E_totalSys)+7),0:18,'Color',rgb('CadetBlue'),'LineStyle','-','LineWidth',3);
hliml = plot((18573.7-workFun/2)*ones(1,length(E_totalSys)+7),0:18,'Color',rgb('CadetBlue'),'LineStyle','-','LineWidth',0.5);
hlimu = plot((18573.7 + workFun/2)*ones(1,length(E_totalSys)+7),0:18,'Color',rgb('CadetBlue'),'LineStyle','-','LineWidth',0.5);
rectangle('Position',[18573.7-workFun/2,0,workFun,15.8],'Facecolor',[rgb('CadetBlue'), 0.5],...
    'EdgeColor','none');
hrec=fill([1 1 0 0],[0 1 1 0],rgb('CadetBlue'));
hold off

axis([18572,18577,0,15.8])
legend([h1,h2,htheo,hrec],{'\color{red} Statistical Uncertainty','Total Uncertainty','Expected Value',...
    'Work Function Uncertainty'},'location','northeast')
legend('boxoff')
PrettyFigureFormat;
yticks(yn)
yticklabels(yana);

%ytickangle(90);
xticks(18572:18576)
xticklabels(["(eV)",e0x]);

title('Effective endpoint including systematics in all analysis types')

%export_fig('plots/E0SummaryAbsolute.pdf','-pdf')