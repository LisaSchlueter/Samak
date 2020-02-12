if ~exist('M','var')
M = MultiRunAnalysis('RunList','KNM1','exclDataStart',14,'chi2','chi2Stat','DataType','Twin','fixPar','5 6 7 8 9 10 11','TwinBias_mnuSq',0);
end
%M.exclDataStart=2;
%M.Fit;
%M.PlotFit('SavePlot','png')

%%
M.exclDataStart=14;
S = RunSensitivity('RunAnaObj',M);
FitResults= S.ComputeFitNuMassTwins('nSamples',1000,'mNuSq',2.2^2);

% M = MultiRunAnalysis('RunList','KNM1','exclDataStart',14,'chi2','chi2Stat','DataType','Twin','fixPar','5 6 7 8 9 10 11','TwinBias_mnuSq',2.2^2);
% M.FitTwin('nSamples',1);
% M.PlotFit('SavePlot','png')
%%
% plot
fig1 = figure(1);
set(fig1, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.7]);
h= histfit(squeeze(FitResults.par(:,1,:)));
PrettyFigureFormat;
h(1).FaceColor = rgb('GoldenRod');
h(1).FaceAlpha = 1;
h(2).Color = rgb('FireBrick'); h(2).LineWidth = 4;
xlabel(sprintf('m_\\nu^2 (eV^2)'));
leg = legend(h(2),sprintf('mean = %.2f eV^2   ,   std = %.2f eV^2',mean(FitResults.par(:,1,:)),std(FitResults.par(:,1,:))));
print(fig1,[getenv('SamakPath'),'knm1ana/knm1Twins/plots/mNuSqFits.png'],'-dpng','-r450');
%%
fig2 = figure(1);
set(fig2, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.7]);
h= histogram(squeeze(FitResults.chi2min)); h.BinWidth=1.75;
hold on;
x = 0:0.1:50;
p = plot(x,chi2pdf(x,FitResults.dof).*numel(FitResults.chi2min).*h.BinWidth);
p.Color = rgb('FireBrick'); p.LineWidth = 4;
PrettyFigureFormat;
h.FaceColor = rgb('SteelBlue');
h.FaceAlpha = 0.6;
hold off;
xlabel(sprintf('\\chi2 (23 dof)'));
print(fig2,[getenv('SamakPath'),'knm1ana/knm1Twins/plots/Chi2Fits.png'],'-dpng','-r450');