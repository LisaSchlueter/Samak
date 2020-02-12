% plot for analysis workshop:
% plot 1 qU distribution and 1 rhod distribution
if ~exist('M','var')
M = MultiRunAnalysis('RunList','KNM1','StackqUCorrFlag','OFF','Debug','ON','DataType','Real');
end

savedir = [getenv('SamakPath'),'knm1ana/knm1_Stacking/plots/'];
MakeDir(savedir);
%%
qu = 16;
f1 = figure('Renderer','opengl');
set(f1, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.6, 0.6]);
h1 = histogram(M.SingleRunData.qU(qu,:)-18575);
h1.FaceColor = rgb('SteelBlue'); h1.FaceAlpha = 0.8;
hold on;
plot((M.RunData.qU(qu,:)-18575).*[1,1],[0,max(h1.BinCounts)+10],'LineWidth',4,'Color',rgb('FireBrick'));
PrettyFigureFormat;
ylim([0,max(h1.BinCounts)+10]);
xlabel('qU - 18575 (eV)')
stdqU = std(M.SingleRunData.qU(qu,:));
leg = legend(sprintf(' single runs \\sigma = %.0f meV',stdqU*1e3),'stacked runs');legend boxoff
print(f1,[savedir,'DataqUDist.png'],'-dpng','-r450');

%%
f2 = figure('Renderer','opengl');
set(f2, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.6, 0.6]);
h2 = histogram(M.SingleRunData.WGTS_MolFrac_TT);
hold on;
plot(M.RunData.WGTS_MolFrac_TT.*[1,1],[0,max(h2.BinCounts)+10],'LineWidth',4,'Color',rgb('FireBrick'));
h2.FaceColor = rgb('GoldenRod'); h2.FaceAlpha = 0.8;
PrettyFigureFormat;
xlabel(sprintf('T_2 mol. frac.'));
ylim([0,max(h2.BinCounts)+10]);
stdTT = std(M.SingleRunData.WGTS_MolFrac_TT);
leg = legend(sprintf('single runs \\sigma = %.1g %%',stdTT/mean(M.SingleRunData.WGTS_MolFrac_TT)),'stacked runs');
legend boxoff
print(f2,[savedir,'DataTTDist.png'],'-dpng','-r450');
