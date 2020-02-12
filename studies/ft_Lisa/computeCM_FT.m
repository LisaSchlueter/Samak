addpath(genpath('../../../Samak2.0'));
myEffects = struct(...
    'RF_EL','ON',...  % Response Function(RF) EnergyLoss
    'RF_BF','ON',...  % RF B-Fields
    'RF_RX','ON',...
    'FSD','ON',...
    'TASR','ON',...
    'TCoff_RAD','ON',...
    'TCoff_OTHER','ON',...
    'DOPoff','ON');  % RF Column Density, Cross Section
    
WGTS_CD_MolPerCm2_RelErr = 0.1;
MACE_Bmax_T_RelErr = 0.02;       % B-Fields
MACE_Ba_T_RelErr= 0.02;
WGTS_B_T_RelErr = 0.02;
ISXsection_RelErr = 0.02;        % Inelastic scattering(IS) cross section  
WGTS_TASR_RelErr=0.023;         % Relative Error on subRun Tritium activity  
FSDNorm_RelErr= 0.03;           % Normalization Error on ground state (sum of GS+ES probability is still always 100%)
FSDShape_RelErr= 0.03;          

%% Model
% RunNr = 40263;
% D = importdata([num2str(RunNr),'.mat']);
% switch RunNr
%     case {40259,40260,40263, 40264, 40265,40266}
%         Data = [D.qU, D.TBDIS, sqrt(D.TBDIS)];
%     case {40257,40258}
%         Data = [D.qU(3:end), D.TBDIS(3:end), sqrt(D.TBDIS(3:end))];
% end
% A = ref_RunSummaries_StackPix(RunNr,'ISCS','Theory','recomputeRF','OFF');  
RunList = [40538:40543,40603,40604,40610:40613,40667:40693];
A = MultiRunAnalysis('RunList',RunList,'chi2','chi2Stat','ringCutFlag','ex2b');
%% Covariance Matrix
WGTS_CD_MolPerCm2_RelErr_local = [0.05, 0.1];
for i=1:numel(WGTS_CD_MolPerCm2_RelErr_local)
    
WGTS_CD_MolPerCm2_RelErr = WGTS_CD_MolPerCm2_RelErr_local(i); 
C = CovarianceMatrix('StudyObject',A.ModelObj,'nTrials',1000,'SysEffect',myEffects,'RecomputeFlag','OFF','SanityPlots','OFF',...                   
                      'WGTS_CD_MolPerCm2_RelErr',WGTS_CD_MolPerCm2_RelErr,...
                      'MACE_Bmax_T_RelErr',MACE_Bmax_T_RelErr,...      
                      'MACE_Ba_T_RelErr',MACE_Ba_T_RelErr,...
                      'WGTS_B_T_RelErr',WGTS_B_T_RelErr,...
                      'ISXsection_RelErr',ISXsection_RelErr,...   
                      'WGTS_TASR_RelErr',WGTS_TASR_RelErr,...
                      'FSDNorm_RelErr',FSDNorm_RelErr,'FSDShape_RelErr',FSDShape_RelErr);
C.ComputeCM_RF; 
%C.ComputeCM; %all ON
end
%CombiCM = C.CovMat;
%CombiCMFrac = C.CovMatFrac;
%C.DecomposeCM('Option','Frac','CombiCM','ON');
%CombiCMFracShape = C.CovMatFracShape;
%CombiCMFracNorm = C.CovMatFracNorm;

%% Plots
% plotOpt ='CM';
% if strcmp(plotOpt,'Frac')
%     plotCM = C.MultiCovMatFrac;
%     plotCombiCM = CombiCMFrac;
%     plotTitle = sprintf('Combined Fractional');
% elseif strcmp(plotOpt,'Shape')
%     plotCM = C.MultiCovMatFracShape;
%     plotCombiCM = CombiCMFracShape;
%     plotTitle = sprintf('Combined Fractional Shape');
% elseif strcmp(plotOpt,'CM')
%     plotCM = C.MultiCovMat;
%     plotCombiCM = CombiCM;
%     plotTitle = sprintf('Combined');
%   elseif strcmp(plotOpt,'Norm')
%     plotCM = C.MultiCovMatFracNorm;
%     plotCombiCM = CombiCMFracNorm;
%     plotTitle = sprintf('Combined Fractional Normalization');  
% end
% %% CM Frac
% fig7  =figure(7);
% set(fig7, 'Units', 'normalized', 'Position', [0.1, 0.1, 1 ,0.8]);
% subplot(2,4,1);
% imagesc(plotCM.CM_RF);
% colorbar;
% PrettyFigureFormat
% title({'Response';' Function'})  
% set(gca,'fontsize', 12);
% pbaspect([1 1 1])
% axis off; 
% 
% subplot(2,4,2);
% imagesc(plotCM.CM_FSD);
% colorbar;
% title({'Final State';' Distribution'});%'FSD')
% PrettyFigureFormat
% set(gca,'fontsize', 12);
% pbaspect([1 1 1])
% axis off; 
% 
% subplot(2,4,5);
% imagesc(plotCM.CM_TCoff);
% colorbar;
% colormap(flipud(gray));
% title({'Theoretical';' Corrections'});
% PrettyFigureFormat
% pbaspect([1 1 1])
% set(gca,'fontsize', 12);
% axis off; 
% 
% subplot(2,4,6);
% imagesc(plotCM.CM_TASR);
% colorbar;
% title({'Activity';' Fluctuation'});
% PrettyFigureFormat
% pbaspect([1 1 1])
% set(gca,'fontsize', 12);
% axis off; 
% 
% subplot(2,4,[3 4 7 8]);
% imagesc(plotCombiCM);
% colorbar;
% colormap(parula);
% title({plotTitle;'Covariance Matrix'});
% PrettyFigureFormat
% pbaspect([1 1 1])
% set(gca,'XTick',[1 8 26],'XTickLabel',{'16.5','18.0','18.6 (kV)'},'FontSize',15,'XMinorTick','off')
% set(gca,'YTick',[3 23],'YTickLabel',{'',''})
% 
% fig100_str = sprintf('./plots/CovarianceMatrix_VFT_%s.eps',plotOpt);
% publish_figure(fig7,fig100_str);
% %% Correlation Plot
% fig8  =figure(8);
%  TickStep = ceil(A.nqU/5);
% set(fig8, 'Units', 'normalized', 'Position', [0.1, 0.1, 1 ,0.9]);
% subplot(2,4,1);
% corplot(plotCM.CM_RF);
% colorbar;
% colormap(flipud(gray));
% PrettyFigureFormat
% %t1 = title({'Response Function';['\Delta \rhod =',num2str(WGTS_CD_MolPerCm2_RelErr),...
%   %  ', \Delta BField =',num2str(MACE_Ba_T_RelErr)];['\Delta \sigma_{inel}=',num2str(ISXsection_RelErr)]},'FontSize',10);
% title({'Response';' Function'})  
% set(gca,'fontsize', 12);
% set(gca,'XTick',1:TickStep:(A.nqU),'XTickLabel',1:TickStep:A.nqU);
%  set(gca,'YTick',1:TickStep:(A.nqU),'YTickLabel',1:TickStep:A.nqU);
% pbaspect([1 1 1])
% 
% subplot(2,4,2);
% corplot(plotCM.CM_FSD);
% colorbar;
% colormap(flipud(gray));
% title('FSD')
% PrettyFigureFormat
% set(gca,'fontsize', 12);
% set(gca,'XTick',1:TickStep:(A.nqU),'XTickLabel',1:TickStep:A.nqU);
%  set(gca,'YTick',1:TickStep:(A.nqU),'YTickLabel',1:TickStep:A.nqU);
% pbaspect([1 1 1])
% 
% subplot(2,4,5);
% corplot(plotCM.CM_TCoff);
% colorbar;
% colormap(flipud(gray));
% title({'Theoretical';' Corrections'});
% PrettyFigureFormat
% pbaspect([1 1 1])
% set(gca,'XTick',1:TickStep:(A.nqU),'XTickLabel',1:TickStep:A.nqU);
% set(gca,'YTick',1:TickStep:(A.nqU),'YTickLabel',1:TickStep:A.nqU);
% set(gca,'fontsize', 12);
% 
% subplot(2,4,6);
% corplot(C.MultiCovMatFrac.CM_TASR);
% colorbar;
% colormap(flipud(gray));
% title({'Activity';' Fluctuation';'per Subrun'});
% PrettyFigureFormat
% pbaspect([1 1 1])
% set(gca,'XTick',1:TickStep:(A.nqU),'XTickLabel',1:TickStep:A.nqU);
% set(gca,'YTick',1:TickStep:(A.nqU),'YTickLabel',1:TickStep:A.nqU);
% set(gca,'fontsize', 12);
% 
% subplot(2,4,[3 4 7 8]);
% corplot(plotCombiCM);
% colorbar;
% colormap(flipud(gray));
% title({'Combined';' Correlation Matrix'});
% PrettyFigureFormat
% pbaspect([1 1 1])
% set(gca,'XTick',1:TickStep:(A.nqU),'XTickLabel',1:TickStep:A.nqU);
% set(gca,'YTick',1:TickStep:(A.nqU),'YTickLabel',1:TickStep:A.nqU);
% %set(gca,'fontsize', 15);
% 
% fig8_str = sprintf('./plots/CorPlot_VFT.pdf');
% publish_figure(fig8,fig8_str);
% 

