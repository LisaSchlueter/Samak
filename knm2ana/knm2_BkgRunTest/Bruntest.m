%
% Run Test Applied to KNM1/ KNM2 scan backgrounds
%

% Retrieve Background Figures
listing = dir(fullfile('../../tritium-data/fit/Knm2/', '*.mat'));

clear B;
for i=1:numel(listing)
    tmp=importdata(listing(i).name);
    B(i)=tmp.B';
    fprintf('File %d = %s - B= %.3f\n',i,listing(i).name,B(i));
end

%% Run Test
Run  = 1:numel(listing);
medB = median(B);
[h,pvalue,stats] = runstest(B,medB,'Alpha',0.05);
%[h,p,stats] = runstest(B,'ud');

%% Plot For Run Test
fig1 = figure('Renderer','opengl');
set(fig1,'units','normalized','pos',[0.1, 0.1,0.9,0.5]);  
p1=plot(Run,B,'s--','Color',rgb('DarkGreen'),'LineWidth',1,'MarkerSize',8,'markerfacecolor','Black');
hold on
l1=line([1 numel(listing)],[medB medB],'Color','Red');
hold off
xlabel('beta-scan');
ylabel('background (cps)');
l=legend([l1,p1],sprintf('median = %.3f cps',medB),sprintf('run test p-value = %.2f (%.0f runs for %.0f variables)',pvalue,stats.nruns,stats.n1+stats.n0));
PrettyFigureFormat

%% Histogram of Scan-Wise Background
fig2 = figure('Renderer','opengl');
set(fig2,'units','normalized','pos',[0.1, 0.1,0.5,0.5]);  
hG=histfit(B);
set(hG(1),'facecolor',rgb('DarkGreen')); set(hG(2),'color','black');set(hG(2),'LineWidth',5)
ylabel('beta-scan');
xlabel('background (cps)');
PrettyFigureFormat

%% Poisson/Gaussian
pdP = fitdist(B','Poisson')
pdG = fitdist(B','Normal')
