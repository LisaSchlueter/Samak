%clear; 
close all;
addpath(genpath('../../../Samak2.0'));

RunList = load('RunList100Good.mat');
RunList = RunList.RunList100Good;
RunList((RunList == 40610) | (RunList == 40612) | (RunList == 40539)) = [];

%if ~exist('SPAll','var')
 %   BuildTable(1);
    load('BuildTableWorkspace')
%end

% All results all runs
ylim = [-2 2];
rr = 2;

%plot(RunList,CD,'ok')
fige0 = figure(8);
pos = get(fige0,'position');
set(fige0,'units','normalized','position',[0.001 0.05 0.85 0.85]);
x = 1:length(RunList);

% PLOT 1
subplot(3,1,1)
E0All_s = SPAll.TableShortStat(:,2) + E0_i - E_SPAll_s;
E0All_s_e = SPAll.TableShortStat(:,6);
hold on
errorb(x,E0All_s,E0All_s_e);
h1 = plot(x,E0All_s,'sk','MarkerFaceColor','black');
plot(x,zeros(1,length(x)),'Color','red','LineStyle','-','LineWidth',3);
hold off
axis([min(x)-1,max(x)+1,ylim]);
xticks(x)
xticklabels([])
ylabel('Short (E_0 - 200 eV)');
title(['<E_0> = ',num2str(round(E_SPAll_s,rr))])
PrettyFigureFormat;

% PLOT 2
subplot(3,1,2)

E0All_m = SPAll.TableMedStat(:,2) + E0_i - E_SPAll_m;
E0All_m_e = SPAll.TableMedStat(:,6);
hold on
errorb(x,E0All_m,E0All_m_e);
h2 = plot(x,E0All_m,'sk','MarkerFaceColor','black');
plot(x,zeros(1,length(x)),'Color','red','LineStyle','-','LineWidth',3);
hold off

axis([min(x)-1,max(x)+1,ylim]);
xticks(x)
xticklabels([]);
ylabel('Medium (E_0 - 400 eV)');
title(['<E_0> = ',num2str(round(E_SPAll_m,rr))])
PrettyFigureFormat;


% PLOT 3
subplot(3,1,3)
E0All_l = SPAll.TableLongStat(:,2) + E0_i - E_SPAll_l;
E0All_l_e = SPAll.TableLongStat(:,6);
hold on
errorb(x,E0All_l,E0All_l_e);
h3 = plot(x,E0All_l,'sk','MarkerFaceColor','black');
plot(x,zeros(1,length(x)),'Color','red','LineStyle','-','LineWidth',3);
hold off
axis([min(x)-1,max(x)+1,ylim]);
ylabel('Long (E_0 - 1600 eV)');
title(['<E_0> = ',num2str(round(E_SPAll_l,rr))])


% ticks for all runs
xticks(x)
xticklabels(string(RunList));
xtickangle(90);
xlabel('Run number');

suptitle('First Tritium - fitted effective E_0 - <E_0> per run (eV)');

PrettyFigureFormat;
export_fig('plots/E0StatAllRuns.pdf','-pdf');
