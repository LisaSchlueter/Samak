addpath(genpath('../../../Samak2.0'));
close all
%BuildTable(1);

 load('BuildTableWorkspace')

%(TotUncVec,StaUncVec,E_totalSys);
yn = 1:12;
E_totalSys = flipud(E_totalSys);
StaUncVec = flipud(StaUncVec);
TotUncVec = flipud(TotUncVec);
 E_totalSys = E_totalSys-E_mean;
StaUncVec = flipud(StaUncVec);
TotUncVec = flipud(TotUncVec);

% yana = ["UAR Short","UAR Med. (*)",...
%     "UAR Long (*)","USR Short (*)",...
%     "USR Med. (*)","USR Long (*)",...
%     "Ring Short","Ring Med. (*)","Single Pix. Short","Single Pix. Med.",...
%     "Multi Pix. Short (*)","Multi Pix. Med. (*)"]';

yana = ["Unseg. Av. Runs Short","Unseg. Av. Runs Med. (*)",...
    "Unseg. Av. Runs Long (*)","Unseg. Stack. Runs Short (*)",...
    "Unseg. Stack. Runs Med. (*)","Unseg. Stack. Runs Long (*)",...
    "Ring Short","Ring Med. (*)","Single Pix. Short","Single Pix. Med.",...
    "Multi Pix. Short (*)","Multi Pix. Med. (*)"]';

e0x = string(-2:2);

yana = flipud(yana);
fige0 = figure(69);
axes('YAxisLocation','right')

pos = get(fige0,'position');
set(fige0,'position',[pos(1:2)/3.5 pos(3:4)*1.5])

hold on
h1 = plot(E_totalSys,yn,'sk','MarkerFaceColor','red');
h2 = plot(E_totalSys,yn,'sk','MarkerFaceColor','black');
errorb(E_totalSys,yn,StaUncVec,'horizontal','color','red');
errorb(E_totalSys,yn,TotUncVec,'horizontal');
htheo = plot(zeros(1,length(E_totalSys)+2),0:13,'Color','blue','LineStyle','--','LineWidth',3);
% hliml = plot(-Unc_mean*ones(1,length(E_totalSys)+3),0:14,'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',3);
% hlimu = plot(Unc_mean*ones(1,length(E_totalSys)+3),0:14,'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',3);


hold off
axis([-2,2,0,14])
% legend([htheo,hliml],{'Theoretical Expectation','Total Uncertainty Unseg. Av. Runs Short Range'},'location','northeast')
legend([h1,h2,htheo],{'\color{red} Statistical Uncertainty','Total uncertainty','Mean Endpoint <E_0>'},'location','northeast')
legend('boxoff')

PrettyFigureFormat;

yticks(yn)
yticklabels(yana);
%ytickangle(90);
xticks(-2:2)
xticklabels(e0x);

title('Relative Endpoint (E_0-<E_0>) including systematics in all analysis types')
export_fig('plots/E0SummaryRelative.pdf','-pdf')