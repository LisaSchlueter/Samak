D = MultiRunAnalysis('RunList','KNM1','FSDFlag','Sibille0p5eV','fitter',...
    'minuit','minuitOpt','min;minos','exclDataStart',14,'chi2','chi2Stat',...
    'TwinBias_mNuSq',0);
S = RunSensitivity('RunAnaObj',D);
%%
nSamples = 5000;
S.ComputeFitNuMassTwins('mNuSq_t',-0.96,'nSamples',nSamples,'ReFit','OFF');
mNuSq = squeeze(S.NFitResults.par(1,1,:));
E0 = squeeze(S.NFitResults.par(1,2,:));

% get rid of outliers
IndexOutlier = find(abs(mNuSq)>30);
mNuSq(IndexOutlier) = [];
E0(IndexOutlier) = [];
%%
sigma = erf(1/sqrt(2));
StatErrLow = prctile(mNuSq,100*(1-sigma)/2)-mean(mNuSq);
StatErrUp  = prctile(mNuSq,(1-(1-sigma)/2)*100)-mean(mNuSq);

StatErrMean = (abs(StatErrLow)+StatErrUp)/2;

%%
f112 = figure('Renderer','opengl');
set(f112, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.8]);

h1 = histogram(mNuSq);
h1.FaceAlpha = 0.9;
h1.BinWidth = 0.18;
h1.FaceColor=rgb('DodgerBlue');
hold on;
h2 = histogram(mNuSq(mNuSq<=StatErrUp+mean(mNuSq) & mNuSq>=StatErrLow+mean(mNuSq)));
h2.BinWidth=h1.BinWidth;h2.FaceColor=rgb('DarkBlue');
PrettyFigureFormat;
xlabel(sprintf('m_\\beta^2 (eV^2)'));
ylabel('experiments');
leg = legend([h1,h2],'5000 pseudo experiments',sprintf('68.3 %%'));
legend box off
set(gca,'FontSize',24);
hold off;

savedir = [getenv('SamakPath'),'/knm1ana/knm1_unblinding/plots/'];
savefile = [savedir,sprintf('FinalFitUncertainty_mNuSqE0_%s_%.0f.png',D.chi2,nSamples)];
print(f112,savefile,'-dpng','-r450');

