close; clear;
index =  [1602  1402  1204  1002   804   602   402   304   202   177   154   129   104    92    82    72    62    52    45    32    22    15     2    -6   -18   -28];
% Stacked Runs Analysis
MRA = MultiRunAnalysis('RunList','StackCD100all','chi2','chi2CM','fixPar','1 5 6','exclDataStart',9);
%MRA.ComputeCM('WGTS_TASR_RelErr',0.001);
MRA.Fit;
E0ErrStackCM = MRA.FitResult.err(2);
MRA_Stat = MultiRunAnalysis('RunList','StackCD100all','chi2','chi2Stat','fixPar','1 5 6','exclDataStart',9);
MRA_Stat.Fit;
E0ErrStackStat = MRA_Stat.FitResult.err(2);

% Runwise Analysis
FitResultsCM   = MRA.FitParameterEvolution;
FitResultsStat = MRA_Stat.FitParameterEvolution;

%%
SigmaSysRunWise = sqrt(FitResultsCM.E0Err.^2-FitResultsStat.E0Err.^2);
SigmaSysStacked = sqrt(E0ErrStackCM^2-E0ErrStackStat^2);
%% Display
figure(1);
s1 = subplot(2,1,1);
x   = linspace(1,numel(MRA.StackedRuns),numel(MRA.StackedRuns));
plot(x,SigmaSysRunWise,'o','color',rgb('CadetBlue'),'MarkerFaceColor',rgb('CadetBlue'))
hold on;
plot(x,repmat(SigmaSysStacked,numel(MRA.StackedRuns),1),'-','LineWidth',2.5,'color',rgb('IndianRed'))
xlabel('run number');
ylabel('\sigma_{sys} on Endpoint (eV)');
leg = legend(['Runwise Fit <\sigma_{sys}> = ',sprintf('%.2f eV, std = %.2f eV',mean(SigmaSysRunWise),std(SigmaSysRunWise))],...
    ['Stacked Runs Fit \sigma_{sys} = ',sprintf('%.2f eV',SigmaSysStacked)],...
    'Location','northwest');
leg.Color = 'none'; legend boxoff;
r = string(MRA.StackedRuns);
xticks(x);
xticklabels(r);
xtickangle(45);
PrettyFigureFormat;
set(gca,'FontSize',18);
set(get(gca,'XAxis'),'FontSize',11);
set(get(gca,'XLabel'),'FontSize',18);
leg.FontSize = 12;
xlim([min(x) max(x)]);
title(sprintf('Samak Endpoint Systematic Uncertainty - %.0f eV below Endpoint',index(MRA.exclDataStart)));

s2 = subplot(2,1,2);
x   = linspace(1,numel(MRA.StackedRuns),numel(MRA.StackedRuns));
yyaxis left 
yAx = get(gca,'YAxis');
set(yAx(1),'Color',rgb('CadetBlue'));
p1 = plot(x,FitResultsCM.E0Err.^2,'o','MarkerFaceColor',rgb('CadetBlue'),'MarkerEdgeColor',rgb('CadetBlue'));
hold on;
ylabel('\sigma_{all}^2 on Endpoint (eV)');
xlabel('run number');
yyaxis right
set(yAx(2),'Color',rgb('CornflowerBlue'));
p2 = plot(x,FitResultsStat.E0Err.^2,'o','MarkerFaceColor',rgb('CornflowerBlue'),'MarkerEdgeColor',rgb('CornflowerBlue'));
ylabel('\sigma_{stat}^2 on Endpoint (eV)');
r = string(MRA.StackedRuns);
xticks(x);
xticklabels(r);
xtickangle(45);
PrettyFigureFormat;
set(gca,'FontSize',18);
set(get(gca,'XAxis'),'FontSize',11);
set(get(gca,'XLabel'),'FontSize',18);
xlim([min(x) max(x)]);

% h = histogram(SigmaSysRunWise,'FaceColor',rgb('CadetBlue'));
% xlabel('\sigma_{sys} on Endpoint (eV)');
% ylabel('number of runs')
% leg = legend(['Runwise Fit <\sigma_{sys}> = ',sprintf('%.2f eV',mean(SigmaSysRunWise))]);
% leg.Color = 'none'; legend boxoff;


