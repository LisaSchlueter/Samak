% plot for KNM1: final fit result
% test FC 
D = MultiRunAnalysis('RunList','KNM1','FSDFlag','Sibille0p5eV','fitter',...
    'minuit','minuitOpt','min;minos','exclDataStart',14,'chi2','chi2Stat',...
    'TwinBias_mNuSq',0);
S = RunSensitivity('RunAnaObj',D);
%%
nSamples = 5000;
S.ComputeFitNuMassTwins('mNuSq_t',0,'nSamples',nSamples,'ReFit','OFF');
mNuSq = S.NFitResults.par(1,1,:);
E0 = S.NFitResults.par(1,2,:);
%mNuSqFit = -0.96;
mNuSqFit = -0.98;
pval = numel(mNuSq(mNuSq<=mNuSqFit))/nSamples*100;
%%
f111 = figure('Renderer','opengl');
set(f111, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.8]);

h1 = histogram(mNuSq);
h1.FaceAlpha = 0.9;
h1.FaceColor=rgb('DodgerBlue');
%h1.EdgeColor = rgb('DarkSlateGray');
h1.BinWidth = 0.245;
hold on;
ymax = max(h1.Values)+20;
h2 = histogram(mNuSq(mNuSq<=mNuSqFit));
h2.FaceColor=rgb('DarkBlue');
h2.BinWidth = h1.BinWidth;
pFit = plot(mNuSqFit.*[1,1],linspace(0,ymax*0.75,2),...%ymax,2),...
    'LineWidth',5,'Color',rgb('Orange'));
ylim([0,ymax]);
PrettyFigureFormat;
xlabel(sprintf('m_\\beta^2 (eV^2)'));
ylabel('experiments');

nu_leg = sprintf('KATRIN best fit: m_{\\beta}^2     =  - %.2g _{-  %.2f}^{+ %.2f} eV^2',...
    0.98,1.06,0.89);%0.96,1.11, 0.92);
pval_leg = sprintf('P ( m_\\beta^2 \\leq - %.2f eV^2 | m_\\beta^2 = 0 eV^2 ) = %.3g %% ',abs(mNuSqFit),pval);
twin_leg = sprintf('pseudo-experiments');
leg = legend([pFit,h1,h2],nu_leg,twin_leg,pval_leg); legend boxoff
leg.Location='northwest';
set(gca,'FontSize',24);
hold off;


savedir = [getenv('SamakPath'),'/knm1ana/knm1_unblinding/plots/'];
savefile = [savedir,sprintf('FinalFitHistogram_mNuSqE0_%s_%.0f.png',D.chi2,nSamples)];
%print(f111,savefile,'-dpng','-r450');
publish_figurePDF(f111,strrep(savefile,'.png','.pdf'));
