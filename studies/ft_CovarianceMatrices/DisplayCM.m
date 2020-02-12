% configuration
TD = 'Flat60';
TimeSec =4*24*60*60;
MACE_Ba_T = 7*1e-04;
MACE_Bmax_T = 6*0.7;
WGTS_B_T    = 3.6*0.7;
WGTS_CD_MolPerCm2 = 5e17;

SysEffect = 'RF';
nTrials = 5000;
FSDNorm_RelErr = 0.01;
FSDShapeGS_RelErr = 0.04;
FSDShapeES_RelErr = 0.18;
MACE_Ba_T_RelErr = 0.02;
MACE_Bmax_T_RelErr = 0.02;
WGTS_B_T_RelErr = 0.02;
WGTS_CD_MolPerCm2_RelErr = 0.05;
ISXsection_RelErr = 0.00;
WGTS_TASR_RelErr = 0.01;

%% Init Model,  compute CM
A= ref_TBD_NominalKATRIN('TD',TD,'TimeSec',TimeSec,'recomputeRF','OFF',...
    'MACE_Ba_T',MACE_Ba_T,'WGTS_B_T',WGTS_B_T,'MACE_Bmax_T',MACE_Bmax_T,...
    'WGTS_CD_MolPerCm2',WGTS_CD_MolPerCm2,'ElossFlag','Abdurashitov');

[CMObj, MultiCM, MultiCMFrac, MultiCMShape, MultiCMNorm] = ...
    ref_CovarianceMatrix_NominalKATRIN('RecomputeFlag','OFF','ModelObj',A,...
    'nTrials',nTrials,'SysEffect',SysEffect,'PlotCM','OFF','SysBudget','99',...
    'FSDNorm_RelErr',FSDNorm_RelErr','FSDShapeGS_RelErr',FSDShapeGS_RelErr,...
    'FSDShapeES_RelErr',FSDShapeES_RelErr,'MACE_Ba_T_RelErr',MACE_Ba_T_RelErr,...
    'MACE_Bmax_T_RelErr',MACE_Bmax_T_RelErr,'WGTS_B_T_RelErr',WGTS_B_T_RelErr,...
    'WGTS_CD_MolPerCm2_RelErr',WGTS_CD_MolPerCm2_RelErr,'ISXsection_RelErr',ISXsection_RelErr,...
    'WGTS_TASR_RelErr',WGTS_TASR_RelErr);
%CMObj.RecomputeFlag = 'ON';
switch SysEffect
    case 'RF'
CMObj.ComputeCM_RF;
    case 'FSD'
        CMObj.ComputeCM_FSD;
    case 'TASR'
        CMObj.ComputeCM_TASR;
    case 'TC'
        CMObj.ComputeCM_TCoff;
    case 'Stacking'
        CMObj.ComputeCM_Stacking;
end
% %% plot
Mode = 'Shape';
switch Mode
    case 'Frac'
        savename = '';
    case 'Shape'
        savename = 'ShapeOnly';
end
%% default Plot
%CMObj.PlotCM('qUWindow',9,'saveplot','ON','PlotEffect',SysEffect,...
 %   'Mode',Mode,'savename',savename,'Convergence','ON');

%% plot correlation matrix
f4 = figure('Renderer','opengl');
set(f4, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7 ,0.8]);
qUWindow = 0;
qUWindowIndex =  min(find(A.qU>=18575-qUWindow));
imagesc(corr(CMObj.CovMatFrac(1:qUWindowIndex,1:qUWindowIndex)));
c = colorbar('northoutside');
colormap(gca,flipud(gray));
c.Label.String =sprintf(['correlation (',SysEffect,')']);
c.FontSize = 22;
PrettyFigureFormat
pbaspect([1 1 1]);
xticks(1:A.nqU);
yticks(1:A.nqU);
grid on
qUmin = sprintf('qU_{min} = E_0- %.0fV',abs(A.qU(1)-18575));
if (A.qU(qUWindowIndex)-18575)>0
    qUmax = sprintf('qU_{max} = E_0+ %.0fV',A.qU(qUWindowIndex)-18575);
elseif  (A.qU(qUWindowIndex)-18575)<0
    qUmax = sprintf('qU_{max} = E_0- %.0fV',abs(A.qU(qUWindowIndex)-18575));
elseif(A.qU(qUWindowIndex)-18575)==0
    qUmax = 'qU_{max} = E_0' ;
end
empty_str1 = strings(10,1);
empty_str2 = strings(A.nqU-23,1);
set(gca,'xticklabel',{empty_str1{:},qUmin,empty_str2{:},qUmax,empty_str1{:}}); set(gca,'yticklabel',[]);
set(gca,'FontSize',22)
            
 f4_str = [sprintf('CorrPlot_ShapeOnlyFrac_%s_%s',TD,SysEffect)];
 publish_figurePDF(f4,['./plots/CovMatInfo/pdf/',f4_str,'.pdf']);
 print(f4,['./plots/CovMatInfo/png/',f4_str,'.png'],'-dpng');
 savefig(f4,['./plots/CovMatInfo/fig/',f4_str,'.fig'],'compact');
 
%% plot cov matrix
f6 = figure('Renderer','opengl');
set(f6, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7 ,0.8]);
        
% Plot Fractional Covariance Matrix
qUWindow = 9;
qUWindowIndex =  min(find(A.qU>=18575-qUWindow));
imagesc(CMObj.CovMatFrac(1:qUWindowIndex,1:qUWindowIndex));
c = colorbar('northoutside');
c.Label.String = sprintf(['frac. covariance (',SysEffect,')']);
c.FontSize = 22;
colormap(parula);
PrettyFigureFormat
pbaspect([1 1 1])
set(gca,'xtick',[7 qUWindowIndex-7]),set(gca,'ytick',[])
qUmin = sprintf('qU_{min} = E_0- %.0fV',abs(A.qU(1)-18575));
if (A.qU(qUWindowIndex)-18575)>0
    qUmax = sprintf('qU_{max} = E_0+ %.0fV',A.qU(qUWindowIndex)-18575);
elseif (A.qU(qUWindowIndex)-18575)<0
    qUmax = sprintf('qU_{max} = E_0- %.0fV',abs(A.qU(qUWindowIndex)-18575));
elseif(A.qU(qUWindowIndex)-18575)==0
    qUmax = 'qU_{max} = E_0' ;
end
set(gca,'xticklabel',{qUmin,qUmax}); set(gca,'yticklabel',[]);
set(gca,'FontSize',22)

 f6_str =[sprintf('CovMat_ShapeOnlyFrac_%s_%s',TD,SysEffect)];
 publish_figurePDF(f6,['./plots/CovMatInfo/pdf/',f6_str,'.pdf']);
 print(f6,['./plots/CovMatInfo/png/',f6_str,'.png'],'-dpng');
 savefig(f6,['./plots/CovMatInfo/fig/',f6_str,'.fig'],'compact');
% %% Plot Convergence Test
% Criteria = 'Cauchy';%'FracVar';%'Determinant';
% [Trials, Convergence] = CMObj.ConvergenceTest('Criterium', Criteria);
% 
% f5 = figure('Renderer','opengl');
% set(f5, 'Units', 'normalized', 'Position', [0.1, 0.1, 1 ,0.5]);
% plot(Trials,Convergence,'Color',rgb('GoldenRod'),'LineWidth',4);
% xlabel('samples');
% switch Criteria
%     case 'Cauchy'
%         ylabel('|| M ||');
%         %xlim([0 5000]);
%     case 'FracVar'
%          ylabel('$\mathbf{\sum_{ij} M_{ij}^{\textrm{frac}}}$','Interpreter','latex');
%          yticks([0.013, 0.015,0.017])
%          %xlim([0 5000]);
%     case 'Determinant'
%          ylabel('det(M)');
% end
% % title(sprintf('Convergence Test'));
% PrettyFigureFormat;
% ylim([min(Convergence) max(Convergence)]);
% set(gca,'FontSize',25);
% grid on;
% 
%  f5_str = [sprintf('Convergence_%s_%s_%s_%.0f',TD,SysEffect,Criteria,nTrials)];
%  publish_figurePDF(f5,['./plots/CovMatInfo/pdf/',f5_str,'.pdf']);
%  print(f5,['./plots/CovMatInfo/png/',f5_str,'.png'],'-dpng');
%  savefig(f5,['./plots/CovMatInfo/fig/',f5_str,'.fig'],'compact');