%TD = 'StackCD100_3hours';
TD = 'StackCD100all';
belowE0 = '402eVbelowE0';
savename = sprintf('./results/AverageRuns_CorrelationCoefficient_RoiPileUp_%sex2_%s.mat',TD,belowE0);
load(savename);
%%
fig13 = figure(13);
set(fig13, 'Units', 'normalized', 'Position', [0.9, 0.9, 0.7, 0.7]);
s1 = subplot(2,1,1);
e1 = errorbar(CorrCoeff,MeanE0,ErrMeanE0,...
    'o','MarkerFaceColor',rgb('CadetBlue'),'Color',rgb('CadetBlue'),'MarkerEdgeColor',rgb('SteelBlue'),'MarkerSize',8,'LineWidth',2);
xlabel('correlation coefficient \rho');
ylabel('<E0_{eff}>');
PrettyFigureFormat;
title(sprintf('weighted mean and standard error of mean \nfor 27 correlated endpoint measurements %s',belowE0));
set(gca,'FontSize',18);

s2 = subplot(2,1,2);
%plot(CorrCoeff, 1-chi2cdf(chi2min,dof(1)), 'o','MarkerFaceColor',rgb('CadetBlue'),'Color',rgb('CadetBlue'),'MarkerSize',8);
plot(CorrCoeff, chi2min, 'o','MarkerFaceColor',rgb('CadetBlue'),'Color',rgb('CadetBlue'),'MarkerEdgeColor',rgb('CornflowerBlue'),'MarkerSize',8);
hold on;
plot(CorrCoeff,dof,'k--','LineWidth',2)
xlabel('correlation coefficient \rho');
%ylabel('p-value (26 dof)');
ylabel(sprintf('\\chi2 (%.0f dof)',dof(1)));
PrettyFigureFormat;
linkaxes([s1,s2],'x');
set(gca,'FontSize',18);
hold off;

savename = sprintf('AverageRuns_CorrelationCoefficient_%s_%s',TD,belowE0);
publish_figurePDF(fig13,['./plots/pdf/',savename,'.pdf']);
export_fig(fig13,['./plots/png/',savename,'.png']);
savefig(fig13,['./plots/fig/',savename,'.fig']);

%%
close;
MRA = MultiRunAnalysis('RunList',TD);
Parameter = 'E0';
Q_i = 18573.7;
PlotPar = ResultsCMall{1}.(Parameter)+Q_i;
PlotUnit = 'eV';

fit_leg = cell(numel(MeanE0),1);
p       = cell(numel(MeanE0)+1,1);
pcolors = colormap(parula);

fig6 = figure(6);
set(fig6, 'Units', 'normalized', 'Position', [0.9, 0.9, 0.9, 0.7]);
x = linspace(1,numel(PlotPar),numel(PlotPar));
p{1} = errorbar(x,PlotPar,ResultsCMall{1}.([Parameter,'Err']),...
    'o','MarkerFaceColor',rgb('CadetBlue'),'Color',rgb('CadetBlue'),'MarkerSize',8);
p{1}.LineWidth = 2;
hold on

for i=1:numel(MeanE0)
 p{i+1} = plot(x,repmat(MeanE0(i),1,numel(x)),...
 '--','LineWidth',2.5,'Color',pcolors((i+5)*3,:));
%fit_leg{i} = sprintf('<%s>= %.2f \\pm %.2f %s\n',...
%    Parameter,MeanE0(i),ErrMeanE0(i),PlotUnit);
fit_leg{i} = ['\rho =',num2str(CorrCoeff(i))];
end
fit_leg{1} = [sprintf('\nMean E0_{eff} \n'),fit_leg{1}];
leg_str = {'Runwise Fit (stat + sys)',fit_leg{:},'location','north'};
leg = legend(leg_str,'Location','bestoutside');
%leg = legend(,fit_leg);
%columnlegend(3,leg_str);
set(leg,'color','none'); legend boxoff;
ylabel(sprintf('%s (%s)',Parameter,PlotUnit));
xlabel('run');
title(sprintf('Samak Fit to %.0f Runs - %s Fit %s - \n Weighted Mean Dependency on Correlation Coefficient',...
    dof(1)+1,Parameter,belowE0));
PrettyFigureFormat;
xlim([1 max(x)]);
ylim([min(PlotPar-ResultsCMall{1}.([Parameter,'Err']))-0.3*max(ResultsCMall{1}.([Parameter,'Err'])) max(PlotPar+ResultsCMall{1}.([Parameter,'Err']))+0.3*max(ResultsCMall{1}.([Parameter,'Err']))]);
set(gca,'FontSize',16);
xticks(x);
xticklabels(string(MRA.StackedRuns));xtickangle(45);
hold off;


savename = sprintf('AverageRuns_MeanE0CorrelationCoefficient_%s_%s',TD,belowE0);
publish_figurePDF(fig6,['./plots/pdf/',savename,'.pdf']);
export_fig(fig6,['./plots/png/',savename,'.png']);
savefig(fig6,['./plots/fig/',savename,'.fig']);

%%
CorrCoeff_idx = find(CorrCoeff==0.9);

fig7 = figure(7);
imagesc(SysCM{CorrCoeff_idx});
colormap('parula');
h = colorbar;
set(gca,'xtick',[1 27]),set(gca,'ytick',[])
set(gca,'xticklabel',{'Run 1',' Run 27'}); set(gca,'yticklabel','');
title(sprintf('Endpoint Covariance Matrix\n  \\rho=%.1f , 27 Runs ',CorrCoeff(CorrCoeff_idx)));
set(get(h,'label'),'string','cov[E0_i, E0_j] (eV²)')
PrettyFigureFormat;
set(gca,'FontSize',18);

savename = sprintf('AverageRuns_SysCME0_%s_%s_CorrCoeff%.1f',TD,belowE0,CorrCoeff(CorrCoeff_idx));
publish_figurePDF(fig7,['./plots/pdf/',savename,'.pdf']);
export_fig(fig7,['./plots/png/',savename,'.png']);
savefig(fig7,['./plots/fig/',savename,'.fig']);

