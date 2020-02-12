clear; close all;

% Add path to samak folder
addpath(genpath('../../../Samak2.0'));

RunList = load('RunListFTFullCD.mat');
RunList = RunList.RunListFTFullCD;
RunList((RunList == 40610) | (RunList == 40612) | (RunList == 40539)) = [];
Nruns = length(RunList);
runs = cell(Nruns,1);
run_name = runs;
qU_M = zeros(26,Nruns);

for r = 1:Nruns
    run_name{r} = [num2str(RunList(r)),'mpix.mat'];
    runs{r} = load(run_name{r});
%     if any(diff(runs{r}.qU(:,1)) < 0)
%        runs{r}.qU(:,1) = sort(runs{r}.qU(:,1));
%     end
    qU_M(:,r) = runs{r}.qU(:,1);
end

%% if runs already loaded run from here
figqU = figure();
qU_mean = mean(qU_M,2);
qU_bias = qU_M - qU_mean;
hold on
h = plot(qU_mean,qU_bias,'o');
set(h, {'MarkerFaceColor'}, get(h,'Color')); 
plot(qU_mean,0.15*ones(1,length(qU_mean)),'Color','red','LineStyle','--','LineWidth',3);
plot(qU_mean,-0.15*ones(1,length(qU_mean)),'Color','red','LineStyle','--','LineWidth',3);
hold off
axis([min(qU_mean)-30,max(qU_mean)+30,-0.23,0.23]);

title('Variation from the mean retarding potential per run');
xlabel('Average retarding potential (V)');
ylabel('Retarding potential variation from the mean (V)');

pos = get(figqU,'position');
set(figqU,'position',[pos(1:2)/3.5 pos(3:4)*1.5])
PrettyFigureFormat;

export_fig('plots/qUStability.pdf','-pdf');




