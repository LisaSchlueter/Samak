% plot for KNM1: correlation plot
D = MultiRunAnalysis('RunList','KNM1','FSDFlag','Sibille0p5eV','fitter',...
    'minuit','minuitOpt','min;minos','exclDataStart',14,'chi2','chi2Stat',...
    'TwinBias_mNuSq',0);
S = RunSensitivity('RunAnaObj',D);
Tweak2Zero = 'OFF'; % shift everything to mNuSq =0
%%
D.chi2 = 'chi2CMShape';
nSamples = 10000;
S.ComputeFitNuMassTwins('mNuSq_t',-0.98,'nSamples',nSamples,'ReFit','OFF');
mNuSq = squeeze(S.NFitResults.par(1,1,:));
E0 = squeeze(S.NFitResults.par(1,2,:));
% get rid of outliers
IndexOutlier = find(abs(mNuSq)>30);
mNuSq(IndexOutlier) = [];
E0(IndexOutlier) = [];
mNuSqFit = -0.98;
E0Fit = 18573.73;
pval = numel(mNuSq(mNuSq<=mNuSqFit))/nSamples*100;

if strcmp(Tweak2Zero,'ON')
   mNuSq = mNuSq+0.98;
   mNuSqFit=0;
end
%% correlation matrix all fit parameters
B = squeeze(S.NFitResults.par(1,3,:));
N = squeeze(S.NFitResults.par(1,4,:));
B(IndexOutlier) = [];
N(IndexOutlier) = [];
CorrMat = corrcoef([mNuSq,E0,B,N]);
cp= corplot(CorrMat);
myticklabels = {sprintf('m^2(\\nu_e)'),sprintf('E_{0,fit}'),'N','B'};
xticklabels(myticklabels);
yticklabels(myticklabels);
PrettyFigureFormat('FontSize',28);
cp.LineWidth= get(gca,'LineWidth');
c = colorbar;
c.Label.String = sprintf('correlation coefficient \\rho');
c.Label.FontSize = get(gca,'FontSize');
c.LineWidth = get(gca,'LineWidth');
c.FontSize = get(gca,'FontSize')-2;

for i=1:size(CorrMat,1)
    for j=1:size(CorrMat,1)
        text(i+0.2,j+0.5,sprintf('%.2f',CorrMat(i,j)),'FontSize',get(gca,'FontSize')-2);
    end
end
set(gca,'XMinorTick','off');
set(gca,'YMinorTick','off');
set(gca,'TickLength',[0 0]);
savedir = [getenv('SamakPath'),'/knm1ana/knm1_unblinding/plots/'];
savefile = [savedir,sprintf('CorrelationMatrix_%s_%.0f.pdf',D.chi2,nSamples)];
publish_figurePDF(gcf,savefile);

%%

group = [zeros(nSamples,1);1];
f111 = figure('Renderer','opengl');
set(f111, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.8]);
p = scatterhist(mNuSq,E0,'Direction','out',...);%,'HistogramDisplayStyle','bar',...
    'Location','Northeast','Color',rgb('DodgerBlue'));
xlabel(sprintf('m^2(\\nu_e) (eV^{ 2})'));
ylabel(sprintf('E_{0,fit} - 18574 (eV)'));
PrettyFigureFormat('FontSize',30);
hold on;
thisYlim = ylim; thisXlim =xlim;
pfit = plot(mNuSqFit,E0Fit-18573.7,'d','MarkerSize',12,'Color',rgb('DarkGoldenRod'),...
    'MarkerFaceColor',rgb('Orange'),'LineWidth',2);

if strcmp(Tweak2Zero,'ON')
    legMode =2;
    str = '0eV2mNuSq';
else
    legMode = 0;
    str='';
end

switch legMode
    case 0
       leg = legend('Pseudo experiments','KATRIN best fit');
    case 1
        % legend with fit results (mNuSq and E0)
        nu_leg = sprintf('m^2(\\nu_e)     =  - %.2f _{-  %.2f}^{+ %.2f} eV^2',...
            round(0.98,0),1.06,0.89);%0.96,1.11, 0.92);
        e0_leg = sprintf('\nE_{0,fit}  =  %.2f _{-  %.2f}^{+ %.2f} eV',...
            18573.73,0.06,0.06);%0.07,0.06);
        leg = legend(pfit,[nu_leg,e0_leg]); legend boxoff;
        leg.Title.String = 'KATRIN best fit';
    case 2
        % simpler legend with correlation coefficient
        leg = legend('Pseudo experiments','Monte Carlo truth');
        %sprintf('correlation coefficient = %.2f',CorrMat(1,2))
end
leg.Location='northwest';
legend boxoff;
savedir = [getenv('SamakPath'),'/knm1ana/knm1_unblinding/plots/'];
savefile = [savedir,sprintf('CorrelationmNuSqE0_%s_%.0f_leg%.0f%s.png',D.chi2,nSamples,legMode,str)];
%export_fig(f111,savefile);

print(f111,savefile,'-dpng','-r100')

