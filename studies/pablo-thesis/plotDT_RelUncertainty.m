clear; close all;

addpath(genpath('../../../Samak2.0'));

% Load RunList
RunList = load('RunList100Good.mat');
RunList = RunList.RunList100Good;
RunList((RunList == 40610) | (RunList == 40612) | (RunList == 40539)) = [];
Nruns = length(RunList);
runs = cell(Nruns,1);
run_name = runs;
DTrelunc = zeros(Nruns,1);
DT = zeros(Nruns,1);
DTabserr = zeros(Nruns,1);
TotalTime = 0;

for r = 1:Nruns
    run_name{r} = [num2str(RunList(r)),'mpix.mat'];
    runs{r} = load(run_name{r});
    DT_subrun = runs{r}.WGTS_MolFrac_DT_SubRun;
    DTrelunc(r) = std(DT_subrun)/mean(DT_subrun);
    DT(r) = runs{r}.WGTS_MolFrac_DT;
    TotalTime = TotalTime + runs{r}.TimeSec;
    DTabserr(r) = runs{r}.WGTS_MolFrac_DT_SubRun_error_mean;
end

DTreluncmean = mean(DTrelunc);
DT_reluncstd = std(DTrelunc);
%plot(RunList,CD,'ok')
figcd = figure(8);
pos = get(figcd,'position');
set(figcd,'position',[pos(1:2)/3.5 pos(3)*2.4 pos(4)*1.12]);
x = 1:length(DTrelunc);

% PLOT 1
subplot(2,1,1)
hold on
errorb(x,DT*100,DTabserr*100);
h1 = plot(x,DT*100,'sk','MarkerFaceColor','black');
%plot(x,0.1*ones(1,length(x)),'Color','red','LineStyle','--','LineWidth',3);
hold off
title('First Tritium - DT concentration (%) ');
axis([min(x)-1,max(x)+1,-inf,inf]);
xticks(x)
xticklabels([])
ylabel('DT concentration (%)');
PrettyFigureFormat;

% PLOT 2
subplot(2,1,2)
hold on
h2 = plot(x,DTrelunc*100,'sk','MarkerFaceColor','black');
%plot(x,0.1*ones(1,length(x)),'Color','red','LineStyle','--','LineWidth',3);
hold off

axis([min(x)-1,max(x)+1,-inf,inf]);
% axis([min(x)-2,max(x)+2,0,0.11]);
curtick = get(gca, 'YTick');
set(gca, 'YTickLabel', cellstr(num2str(curtick(:))));
xticks(x)
xticklabels(string(RunList));
xtickangle(90);

title('First Tritium - DT relative uncertainty in each run');
ylabel(sprintf('Relative uncertainty %% \n(\\Delta DT/DT)*100'),'interp','tex');

xlabel('Run number');


PrettyFigureFormat;
export_fig('plots/DTRelUncSubRun.pdf','-pdf');


