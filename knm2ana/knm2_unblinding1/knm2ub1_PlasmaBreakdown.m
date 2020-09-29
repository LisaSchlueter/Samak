% Calculte and/or plot systematics breakdown for KNM2 unblinding 1


file  = [getenv('SamakPath'),'inputs/CovMat/LongPlasma/CM/LongPlasma_KNM2Prompt_5000Trials_0meVEloss_99000meVElossErr_122meVSigma_127meVSigmaErr_0CorrCoeff_Troitsk+.mat'];
d = importdata(file);

%%
plotpath = [getenv('SamakPath'),'knm2ana/knm2_unblinding1/plots/'];
MakeDir(plotpath)
plotnameBroadening = [plotpath,'knm2ub1_PlasmaBreakdown_Broadening.png'];

f1 = figure('Units','normalized','Position',[0.1,0.1,0.4,0.45]);
h1 = histogram(d.MACE_Var_v,'Normalization','probability','FaceColor',rgb('SkyBlue'),'FaceAlpha',1);
PrettyFigureFormat('FontSize',22);
xlabel(sprintf('\\sigma_0^2 (eV^2)'));
ylabel('Frequency');
set(gca,'LineWidth',1);
legend(sprintf('\\mu = %.4f eV^{ 2} , \\sigma = %.4f eV^{ 2}',mean(d.MACE_Var_v),std(d.MACE_Var_v)),'EdgeColor',rgb('Silver'),'Location','northwest');
xlim([-0.05 0.08])
t = title(sprintf('%.0f samples',numel(d.is_EOffset_v)));
t.FontWeight = 'normal'; t.FontSize = get(gca,'FontSize');

print(f1,plotnameBroadening,'-dpng','-r350');
%%
f2 = figure('Units','normalized','Position',[0.1,0.1,0.4,0.45]);
h1 = histogram(d.is_EOffset_v,'Normalization','probability','FaceColor',rgb('GreenYellow'),'FaceAlpha',1);
PrettyFigureFormat('FontSize',22);
xlabel(sprintf('\\Delta_{10} (eV)'));
ylabel('Frequency');
set(gca,'LineWidth',1);
legend(sprintf('\\mu = %.4f eV , \\sigma = %.4f eV',mean(d.is_EOffset_v),std(d.is_EOffset_v)),'EdgeColor',rgb('Silver'),'Location','northwest');
t = title(sprintf('%.0f samples',numel(d.is_EOffset_v)));
t.FontWeight = 'normal'; t.FontSize = get(gca,'FontSize');

plotnameShift = [plotpath,'knm2ub1_PlasmaBreakdown_ElossShift.png'];
print(f2,plotnameShift,'-dpng','-r350');

