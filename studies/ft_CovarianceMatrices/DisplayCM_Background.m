TD = 'StackCD100allex2';
TimeSec =45*24*60*60*(124/148);
MACE_Ba_T = 7*1e-04;
MACE_Bmax_T = 6*0.7;
WGTS_B_T    = 3.6*0.7;
WGTS_CD_MolPerCm2 = 5e17;
BKG_RateSec = 0.350;
A= ref_TBD_NominalKATRIN('TD',TD,'TimeSec',TimeSec,'recomputeRF','OFF',...
    'MACE_Ba_T',MACE_Ba_T,'WGTS_B_T',WGTS_B_T,'MACE_Bmax_T',MACE_Bmax_T,...
    'WGTS_CD_MolPerCm2',WGTS_CD_MolPerCm2,'BKG_RateAllFPDSec',BKG_RateSec);
%%

[CovMat, CovMatFrac, CovMatShape, CovMatFracShape] = ComputeCM_Background(...
'StudyObject',A,'plotFit','ON','RecomputeFlag','ON','nTrials',5500,'savePlot','ON');

%% corr plot
f4 = figure('Renderer','opengl');
set(f4, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7 ,0.8]);
imagesc(corr(CovMatFracShape));
c = colorbar('northoutside');
colormap(gca,flipud(gray));
c.Label.String =sprintf(['correlation (bkg)']);
c.FontSize = 22;
PrettyFigureFormat
pbaspect([1 1 1])
xticks(1:A.nqU);
yticks(1:A.nqU);
grid on

qUmin = sprintf('qU_{min} = E_0- %.0fV',abs(A.qU(1)-18575));
qUmax = sprintf('qU_{max} = E_0+ %.0fV',abs(A.qU(end)-18575));
empty_str1 = strings(4,1);
empty_str2 = strings(A.nqU-10,1);
set(gca,'xticklabel',{empty_str1{:},qUmin,empty_str2{:},qUmax,empty_str1{:}}); set(gca,'yticklabel',[]);
set(gca,'FontSize',22)
            
 f4_str = [sprintf('CorrPlot_ShapeOnlyFrac_%s_BKG',TD)];
 publish_figurePDF(f4,['./plots/CovMatInfo/pdf/',f4_str,'.pdf']);
 print(f4,['./plots/CovMatInfo/png/',f4_str,'.png'],'-dpng');
 savefig(f4,['./plots/CovMatInfo/fig/',f4_str,'.fig'],'compact');
 
%%
f6 = figure('Renderer','opengl');
set(f6, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7 ,0.8]);
        
% Plot Fractional Covariance Matrix
imagesc(CovMatFracShape);
c = colorbar('northoutside');
c.Label.String = sprintf(['frac. covariance (bkg)']);
c.FontSize = 22;
colormap(parula);
PrettyFigureFormat
pbaspect([1 1 1])
set(gca,'xtick',[1,A.nqU]),set(gca,'ytick',[])
qUmin = sprintf('qU_{min} = E_0- %.0fV',abs(A.qU(1)-18575));
qUmax = sprintf('qU_{max} = E_0+ %.0fV',A.qU(end)-18575);
set(gca,'xticklabel',{qUmin,qUmax}); set(gca,'yticklabel',[]);
set(gca,'FontSize',22)

 f6_str =[sprintf('CovMat_ShapeOnlyFrac_%s_Bkg',TD)];
 publish_figurePDF(f6,['./plots/CovMatInfo/pdf/',f6_str,'.pdf']);
 print(f6,['./plots/CovMatInfo/png/',f6_str,'.png'],'-dpng');
 savefig(f6,['./plots/CovMatInfo/fig/',f6_str,'.fig'],'compact');
 
 %%