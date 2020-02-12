% d = load('CMStack-35Runs_RunList-StackCD100downex2_1000Trials.mat');
% CM = d.obj.CovMat;
% CMfrac = d.obj.CovMatFrac;


%%


M = MultiRunAnalysis('RunList','StackCD100all','chi2','chi2Stat');
nTrials  = 5000;
M.ComputeCM('nTrials',nTrials);

%M.StackCM_Obj.PlotCM;
%% plot correlation matrix
f4 = figure('Renderer','opengl');
set(f4, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7 ,0.8]);
qUWindow = 0;
qUWindowIndex =  min(find(M.ModelObj.qU>=18575-qUWindow));
imagesc(corr(M.StackCM_Obj.CovMatFrac(1:qUWindowIndex,1:qUWindowIndex)));
c = colorbar('northoutside');
colormap(gca,flipud(gray));
c.Label.String ='correlation (Stacking)';
c.FontSize = 22;
PrettyFigureFormat
pbaspect([1 1 1])
xticks(1:M.ModelObj.nqU);
yticks(1:M.ModelObj.nqU);
grid on
qUmin = sprintf('qU_{min} = E_0- %.0fV',abs(M.ModelObj.qU(1)-18575));
qUmax = 'qU_{max} = E_0';
empty_str1 = strings(4,1);
empty_str2 = strings(M.ModelObj.nqU-10,1);
set(gca,'xticklabel',{empty_str1{:},qUmin,empty_str2{:},qUmax,empty_str1{:}}); set(gca,'yticklabel',[]);
set(gca,'FontSize',22)

f4_str = [sprintf('CorrPlot_ShapeOnlyFrac_%s_Stacking',M.ModelObj.TD)];
publish_figurePDF(f4,['./plots/CovMatInfo/pdf/',f4_str,'.pdf']);
print(f4,['./plots/CovMatInfo/png/',f4_str,'.png'],'-dpng');
savefig(f4,['./plots/CovMatInfo/fig/',f4_str,'.fig'],'compact');
%% plot cov matrix
f6 = figure('Renderer','opengl');
set(f6, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7 ,0.8]);
        
% Plot Fractional Covariance Matrix
qUWindow = 25;
qUWindowIndex =  min(find(M.ModelObj.qU>=18575-qUWindow));
imagesc(M.StackCM_Obj.CovMatFrac(1:qUWindowIndex,1:qUWindowIndex));
c = colorbar('northoutside');
c.Label.String = 'frac. covariance (Stacking)';
c.FontSize = 22;
colormap(parula);
PrettyFigureFormat
pbaspect([1 1 1])
set(gca,'xtick',[5 qUWindowIndex-3]),set(gca,'ytick',[])
qUmin = sprintf('qU_{min} = E_0- %.0fV',abs(M.ModelObj.qU(1)-18575));
if (M.ModelObj.qU(qUWindowIndex)-18575)>0
    qUmax = sprintf('qU_{max} = E_0+ %.0fV',M.ModelObj.qU(qUWindowIndex)-18575);
elseif (M.ModelObj.qU(qUWindowIndex)-18575)<0
    qUmax = sprintf('qU_{max} = E_0- %.0fV',abs(M.ModelObj.qU(qUWindowIndex)-18575));
elseif(M.ModelObj.qU(qUWindowIndex)-18575)==0
    qUmax = 'qU_{max} = E_0' ;
end
set(gca,'xticklabel',{qUmin,qUmax}); set(gca,'yticklabel',[]);
set(gca,'FontSize',22)

 f6_str =[sprintf('CovMat_ShapeOnlyFrac_%s_%s',M.ModelObj.TD,'Stacking')];
 publish_figurePDF(f6,['./plots/CovMatInfo/pdf/',f6_str,'.pdf']);
 print(f6,['./plots/CovMatInfo/png/',f6_str,'.png'],'-dpng');
 savefig(f6,['./plots/CovMatInfo/fig/',f6_str,'.fig'],'compact');


 %%
 Criteria = 'Cauchy';%'FracVar';%'Determinant';
[Trials, Convergence] = M.StackCM_Obj.ConvergenceTest('Criterium', Criteria);

f5 = figure('Renderer','opengl');
set(f5, 'Units', 'normalized', 'Position', [0.1, 0.1, 1 ,0.5]);
plot(Trials,Convergence,'Color',rgb('GoldenRod'),'LineWidth',4);
xlabel('samples');
switch Criteria
    case 'Cauchy'
        ylabel('|| M ||');
        %xlim([0 5000]);
    case 'FracVar'
         ylabel('$\mathbf{\sum_{ij} M_{ij}^{\textrm{frac}}}$','Interpreter','latex');
         yticks([0.013, 0.015,0.017])
         %xlim([0 5000]);
    case 'Determinant'
         ylabel('det(M)');
end
PrettyFigureFormat;
ylim([min(Convergence) max(Convergence)]);
set(gca,'FontSize',25);
grid on;

 f5_str = [sprintf('Convergence_%s_Stacking_%s_%.0f',M.ModelObj.TD,Criteria,nTrials)];
 publish_figurePDF(f5,['./plots/CovMatInfo/pdf/',f5_str,'.pdf']);
 print(f5,['./plots/CovMatInfo/png/',f5_str,'.png'],'-dpng');
 savefig(f5,['./plots/CovMatInfo/fig/',f5_str,'.fig'],'compact');