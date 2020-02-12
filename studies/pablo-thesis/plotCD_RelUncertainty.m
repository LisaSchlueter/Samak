clear; close all;

addpath(genpath('../../../Samak2.0'));

% Load RunList
RunList = load('RunList100Good.mat');
RunList = RunList.RunList100Good;
RunList((RunList == 40610) | (RunList == 40612) | (RunList == 40539)) = [];
Nruns = length(RunList);
runs = cell(Nruns,1);
run_name = runs;
CDrelerr = zeros(Nruns,1);
CD = zeros(Nruns,1);
TotalTime = 0;

for r = 1:Nruns
    run_name{r} = [num2str(RunList(r)),'mpix.mat'];
    runs{r} = load(run_name{r});
    CD(r) = runs{r}.WGTS_CD_MolPerCm2;
    CD_subrun = runs{r}.WGTS_CD_MolPerCm2_SubRun;
    CDrelerr(r) = std(CD_subrun)/mean(CD_subrun);
    TotalTime = TotalTime + runs{r}.TimeSec;
end

CDreluncmean = mean(CDrelerr);
CD_reluncstd = std(CDrelerr);

%plot(RunList,CD,'ok')
figcd = figure(8);
pos = get(figcd,'position');
set(figcd,'position',[pos(1:2)/3.5 pos(3)*2.4 pos(4)*1.12]);
x = 1:length(CDrelerr);

% PLOT 1
subplot(2,1,1)
hold on
errorb(x,CD,CDrelerr);
h1 = plot(x,CD,'sk','MarkerFaceColor','black');
%plot(x,0.1*ones(1,length(x)),'Color','red','LineStyle','--','LineWidth',3);
hold off
title('First Tritium - Column Density (molecules/cm^2)');
axis([min(x)-1,max(x)+1,-inf,inf]);
xticks(x)
xticklabels([])
ylabel(sprintf('Column Density \n(molecules/cm^2)'));
PrettyFigureFormat;

% PLOT 2
subplot(2,1,2)
hold on
h2 = plot(x,CDrelerr,'sk','MarkerFaceColor','black');
%plot(x,0.1*ones(1,length(x)),'Color','red','LineStyle','--','LineWidth',3);
hold off

axis([min(x)-1,max(x)+1,-inf,inf]);
% axis([min(x)-2,max(x)+2,0,0.11]);
curtick = get(gca, 'YTick');
set(gca, 'YTickLabel', cellstr(num2str(curtick(:))));
xticks(x)
xticklabels(string(RunList));
xtickangle(90);

title('First Tritium - Column Density relative uncertainty ');
ylabel(sprintf('Relative uncertainty %% \n(\\Delta\\rho_d/\\rho_d)*100'),'interp','tex');
xlabel('Run number');


PrettyFigureFormat;
export_fig('plots/CDRelUncSubRun.pdf','-pdf');


