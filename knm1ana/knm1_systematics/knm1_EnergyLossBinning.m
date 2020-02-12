%% compute covariance matices for sys. uncertainty on energy loss function
% use different energy loss function binnings (different energy intervals)
% compare covariance matrices: residuals and impact on nu-mass sensitivity

ELossFlag = 'KatrinT2';%'Abdurashitov';%''Aseev';%Abdurashitov';%KatrinD2';%'Abdurashitov';%'KatrinD2';
nTrials = 1000;
RecomputeFlag = 'OFF';

% Init Model Object and covariance matrix object
if ~exist('A','var')
    A = MultiRunAnalysis('RunList','KNM1','chi2','chi2Stat','DataType','Twin');
end


ELossBinStep = 0.2;
RFBinStep = 0.04;

A.ModelObj.ELossFlag = ELossFlag; A.ModelObj.InitializeELossFunction;
SysEffects = struct('RF_EL','ON');
CM = CovarianceMatrix('StudyObject',A.ModelObj, 'nTrials',nTrials,'SysEffect',SysEffects,'RecomputeFlag',RecomputeFlag);

maxE_el = 3000;
CM.ComputeCM_RF('maxE_el',maxE_el,'ELossBinStep',ELossBinStep,'RFBinStep',RFBinStep);
CovMatFrac3000 = CM.MultiCovMatFrac.CM_RF_EL;
CovMat3000 = CM.MultiCovMat.CM_RF_EL;

maxE_el = 2000;
CM.ComputeCM_RF('maxE_el',maxE_el,'ELossBinStep',ELossBinStep,'RFBinStep',RFBinStep);
CovMatFrac2000 = CM.MultiCovMatFrac.CM_RF_EL;
CovMat2000 = CM.MultiCovMat.CM_RF_EL;

maxE_el = 1000;
CM.ComputeCM_RF('maxE_el',maxE_el,'ELossBinStep',ELossBinStep,'RFBinStep',RFBinStep);
CovMatFrac1000 = CM.MultiCovMatFrac.CM_RF_EL;
CovMat1000 = CM.MultiCovMat.CM_RF_EL;

maxE_el = 500;
CM.ComputeCM_RF('maxE_el',maxE_el,'ELossBinStep',ELossBinStep,'RFBinStep',RFBinStep);
CovMatFrac500 = CM.MultiCovMatFrac.CM_RF_EL;
CovMat500 = CM.MultiCovMat.CM_RF_EL;

maxE_el = 200;
CM.ComputeCM_RF('maxE_el',maxE_el,'ELossBinStep',ELossBinStep,'RFBinStep',RFBinStep);
CovMatFrac200 = CM.MultiCovMatFrac.CM_RF_EL;
CovMat200 = CM.MultiCovMat.CM_RF_EL;
%%
f33 = figure('Renderer','opengl');
set(f33, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7,1]);

subplot(2,2,1);
imagesc(CovMatFrac3000(2:end,2:end)-CovMatFrac2000(2:end,2:end))
title('residuals: 3keV - 2keV')
PrettyFigureFormat;
c = colorbar;
c.Label.String = 'fractional covariance';
c.FontSize = 20;
colormap(parula);
PrettyFigureFormat
pbaspect([1 1 1])
set(gca,'xtick',[3 37]),set(gca,'ytick',[])
set(gca,'xticklabel',{[num2str(round(A.ModelObj.qU(2))),'eV'],[num2str(round(A.ModelObj.qU(end))),'eV']}); set(gca,'yticklabel',[]);

subplot(2,2,2);
imagesc(CovMatFrac3000(2:end,2:end)-CovMatFrac1000(2:end,2:end))
title('residuals: 3keV - 1keV')
PrettyFigureFormat;
c = colorbar;
c.Label.String = 'fractional covariance';
c.FontSize = 20;
colormap(parula);
PrettyFigureFormat
pbaspect([1 1 1])
set(gca,'xtick',[3 37]),set(gca,'ytick',[])
set(gca,'xticklabel',{[num2str(round(A.ModelObj.qU(2))),'eV'],[num2str(round(A.ModelObj.qU(end))),'eV']}); set(gca,'yticklabel',[]);

subplot(2,2,3);
imagesc(CovMatFrac3000(2:end,2:end)-CovMatFrac500(2:end,2:end))
title('residuals: 3keV - 0.5keV')
PrettyFigureFormat;
c = colorbar;
c.Label.String = 'fractional covariance';
c.FontSize = 20;
colormap(parula);
PrettyFigureFormat
pbaspect([1 1 1])
set(gca,'xtick',[3 37]),set(gca,'ytick',[])
set(gca,'xticklabel',{[num2str(round(A.ModelObj.qU(2))),'eV'],[num2str(round(A.ModelObj.qU(end))),'eV']}); set(gca,'yticklabel',[]);

subplot(2,2,4);
imagesc(CovMatFrac3000(2:end,2:end)-CovMatFrac200(2:end,2:end))
title('residuals: 3keV - 0.2keV')
PrettyFigureFormat;
c = colorbar;
c.Label.String = 'fractional covariance';
c.FontSize = 20;
colormap(parula);
PrettyFigureFormat
pbaspect([1 1 1])
set(gca,'xtick',[3 37]),set(gca,'ytick',[])
set(gca,'xticklabel',{[num2str(round(A.ModelObj.qU(2))),'eV'],[num2str(round(A.ModelObj.qU(end))),'eV']}); set(gca,'yticklabel',[]);

% save plot
savepath = [getenv('SamakPath'),'knm1ana/knm1_systematics/plots/'];
savename = [savepath,sprintf('knm1_eloss_%s_CM_Binning.png',A.ModelObj.ELossFlag)];

if ~exist(savepath,'dir')
    system(['mkdir ',savepath]);
end
print(gcf,savename,'-dpng','-r400');




%%



