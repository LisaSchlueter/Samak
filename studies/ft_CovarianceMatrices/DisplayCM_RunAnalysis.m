R = MultiRunAnalysis('RunList','StackCD100all','chi2','chi2Stat','DataEffCor','RunSummary');
R.ComputeCM('nTrials',5000,'Stack','ON');
%%
R.FitCM_Obj.SysEffect = struct('RF_EL','ON','RF_BF','ON','RF_RX','ON');
R.FitCM_Obj.ComputeCM_RF;
SysEffect = 'RF';
f6 = figure('Renderer','opengl');
set(f6, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7 ,0.8]);
        
% Plot Fractional Covariance Matrix
qUWindow = 2;
qUWindowIndex =  26;%min(find(R.ModelObj.qU>=18575-qUWindow));
imagesc(R.FitCM_Obj.CovMatFracShape(1:qUWindowIndex,1:qUWindowIndex));
c = colorbar('northoutside');
c.Label.String = sprintf(['frac. covariance (',SysEffect,')']);
c.FontSize = 22;
colormap(parula);
PrettyFigureFormat
pbaspect([1 1 1])
set(gca,'xtick',[6 qUWindowIndex-4]),set(gca,'ytick',[])
qUmin = sprintf('qU_{min} = E_0- %.0fV',abs(R.ModelObj.qU(1)-18575));
if (R.ModelObj.qU(qUWindowIndex)-18575)>=0
    qUmax = sprintf('qU_{max} = E_0+ %.0fV',R.ModelObj.qU(qUWindowIndex)-18575);
else
    qUmax = sprintf('qU_{max} = E_0- %.0fV',abs(R.ModelObj.qU(qUWindowIndex)-18575));
end
set(gca,'xticklabel',{qUmin,qUmax}); set(gca,'yticklabel',[]);
set(gca,'FontSize',22)

 f6_str =[sprintf('CovMat_ShapeOnlyFrac_%s_%s',R.ModelObj.TD,SysEffect)];
 publish_figurePDF(f6,['./plots/CovMatInfo/pdf/',f6_str,'.pdf']);
 print(f6,['./plots/CovMatInfo/png/',f6_str,'.png'],'-dpng');
 savefig(f6,['./plots/CovMatInfo/fig/',f6_str,'.fig'],'compact');