%% KNM-1 Model
chi2          = 'chi2CMShape';
DataType      = 'Real';
KNM1SysBudget = 24;
KNM1Doppler   = 'OFF';
savedir = [getenv('SamakPath'),'knm2ana/knm2_Combination/results/'];
savename = sprintf('%sknm2_CombiFit_%s_Uniform_%s_Knm1SysBudget%.0f_Knm1DE%s.mat',savedir,DataType,chi2,KNM1SysBudget,KNM1Doppler);

if exist(savename,'file') 
    d = importdata(savename);
    fprintf('load results from file %s \n',savename)
else
    fprintf(2, 'file not found %s \n',savename);
    return
end

%% prepare cov mat

nqU1 = 27;
nqU2 = 28;
nqU = nqU1+nqU2;
COVMAT          = diag(zeros(nqU,nqU));
COVMATShape     = diag(zeros(nqU,nqU));
COVMATFracShape = diag(zeros(nqU,nqU));

COVMAT(1:nqU1,1:nqU1) = d.K1.FitCM(d.K1.exclDataStart:end,d.K1.exclDataStart:end);
COVMAT(nqU1+1:nqU,nqU1+1:nqU) = d.K2.FitCM(d.K2.exclDataStart:end,d.K2.exclDataStart:end);

COVMATShape(1:nqU1,1:nqU1) = d.K1.FitCMShape(d.K1.exclDataStart:end,d.K1.exclDataStart:end);
COVMATShape(nqU1+1:nqU,nqU1+1:nqU) = d.K2.FitCMShape(d.K2.exclDataStart:end,d.K2.exclDataStart:end);

COVMATFracShape(1:nqU1,1:nqU1) = d.K1.FitCMFracShape(d.K1.exclDataStart:end,d.K1.exclDataStart:end);
COVMATFracShape(nqU1+1:nqU,nqU1+1:nqU) = d.K2.FitCMFracShape(d.K2.exclDataStart:end,d.K2.exclDataStart:end);

%% cov mat shape-only
close all
imagesc(COVMATShape)
cb = colorbar;
PrettyFigureFormat('FontSize',20)
pbaspect([1 1 1]);
cb.Label.String = sprintf('Covariance (shape-only)');
cb.Label.FontSize = get(gca,'FontSize')+2;
xlabel('Retarding potential bin');
ylabel('Retarding potential bin');
plotdir = [getenv('SamakPath'),'knm2ana/knm2_Combination/plots/'];
plotname = sprintf('%sknm2_CombiFitUniform_CovMatShape.png',plotdir);
print(gcf,plotname,'-dpng','-r300');

%% cov mat shape-only zoom
close all
imagesc(COVMATShape(1:nqU1,1:nqU1));
cb = colorbar;
PrettyFigureFormat('FontSize',20)
pbaspect([1 1 1]);
cb.Label.String = sprintf('Covariance (shape-only)');
cb.Label.FontSize = get(gca,'FontSize')+2;
xlabel('Retarding potential bin');
ylabel('Retarding potential bin');
plotdir = [getenv('SamakPath'),'knm2ana/knm2_Combination/plots/'];
plotname = sprintf('%sknm2_CombiFitUniform_CovMatShapeKNM1Zoom.png',plotdir);
print(gcf,plotname,'-dpng','-r300');

%% cov mat frac shape-only
close all
imagesc(COVMATFracShape)
cb = colorbar;
PrettyFigureFormat('FontSize',20)
pbaspect([1 1 1]);
cb.Label.String = sprintf('Covariance (shape-only)');
cb.Label.FontSize = get(gca,'FontSize')+2;
xlabel('Retarding potential bin');
ylabel('Retarding potential bin');
plotdir = [getenv('SamakPath'),'knm2ana/knm2_Combination/plots/'];
plotname = sprintf('%sknm2_CombiFitUniform_CovMatFracShape.png',plotdir);
print(gcf,plotname,'-dpng','-r300');

%% cov mat shape-only zoom
close all
imagesc(COVMATFracShape(1:nqU1,1:nqU1));
cb = colorbar;
PrettyFigureFormat('FontSize',20)
pbaspect([1 1 1]);
cb.Label.String = sprintf('Covariance (shape-only)');
cb.Label.FontSize = get(gca,'FontSize')+2;
xlabel('Retarding potential bin');
ylabel('Retarding potential bin');
plotdir = [getenv('SamakPath'),'knm2ana/knm2_Combination/plots/'];
plotname = sprintf('%sknm2_CombiFitUniform_CovMatFracShapeKNM1Zoom.png',plotdir);
print(gcf,plotname,'-dpng','-r300');