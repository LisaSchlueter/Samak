
%Display several Covariance Matrices
MRA = MultiRunAnalysis('RunList','StackCD100all','chi2','chi2Stat');

WGTS_TASR_RelErr = 0.001;
FSDNorm_RelErr = 0.01;
FSDShapeGS_RelErr = 0.04;
FSDShapeES_RelErr = 0.18;
WGTS_CD_MolPerCm2_RelErr = 0.08;
nTrials = 1000;
saveplot = 'ON';

MRA.ComputeCM('nTrials',nTrials,'RecomputeFlag','OFF',...
    'WGTS_TASR_RelErr',WGTS_TASR_RelErr,...
    'WGTS_CD_MolPerCm2_RelErr',WGTS_CD_MolPerCm2_RelErr,...
    'FSDNorm_RelErr',FSDNorm_RelErr,'FSDShapeGS_RelErr',FSDShapeGS_RelErr,'FSDShapeES_RelErr',FSDShapeES_RelErr);
%% plot
f5 = figure(5);
set(f5, 'Units', 'normalized', 'Position', [0.1, 0.1, 1, 0.8]);
subplot(2,4,1)
imagesc(MRA.FitCM_Obj.MultiCovMatFrac.CM_RF); pbaspect([1 1 1]);
colorbar;
PrettyFigureFormat;
title(sprintf('Response \nFunction'),'FontSize',26);
set(gca,'xtick',[]),set(gca,'ytick',[])
set(gca,'xticklabel',[]), set(gca,'yticklabel',[]);

subplot(2,4,2)
imagesc(MRA.FitCM_Obj.MultiCovMatFrac.CM_FSD); pbaspect([1 1 1]);
colorbar;
PrettyFigureFormat;
title(sprintf('Final State \nDistribution'),'FontSize',26);
set(gca,'xtick',[]),set(gca,'ytick',[])
set(gca,'xticklabel',[]), set(gca,'yticklabel',[]);

subplot(2,4,5)
imagesc(MRA.FitCM_Obj.MultiCovMatFrac.CM_TASR); pbaspect([1 1 1]);
colorbar;
PrettyFigureFormat;
title(sprintf('Tritium Activity \nFluctuation'),'FontSize',26);
set(gca,'xtick',[]),set(gca,'ytick',[])
set(gca,'xticklabel',[]), set(gca,'yticklabel',[]);

subplot(2,4,6)
imagesc(MRA.FitCM_Obj.MultiCovMatFrac.CM_TCoff); pbaspect([1 1 1]);
colorbar;
PrettyFigureFormat;
title(sprintf('Theoretical \nCorrections'),'FontSize',26);
set(gca,'xtick',[]),set(gca,'ytick',[])
set(gca,'xticklabel',[]), set(gca,'yticklabel',[]);

subplot(2,4,[3,4,7,8])
imagesc(MRA.FitCM_Obj.CovMatFrac); pbaspect([1 1 1]);
colorbar;
PrettyFigureFormat;
title(sprintf('Combined Fractional \nCovariance Matrix'),'FontSize',26);
set(gca,'xtick',[1 length(MRA.FitCM_Obj.CovMatFrac)]),set(gca,'ytick',[])
set(gca,'xticklabel',{'qU_{min}',' qU_{max}'}); set(gca,'yticklabel',[]);
set(gca,'FontSize',26);

if strcmp(saveplot,'ON')
    if  exist('plots/')==0 %if folder doens't exist
        mkdir plots
    end
    save_name = sprintf('./plots/CovMatOverview_%s',MRA.ModelObj.TD);
    publish_figurePDF(f5,[save_name,'.pdf']);
    export_fig(f5,[save_name,'.png']);
    savefig(f5,[save_name,'.fig']);
end