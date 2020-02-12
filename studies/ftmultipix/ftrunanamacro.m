%clear; close all;
addpath(genpath('../../../Samak2.0'));
%40538,40539,40540,40541,40542,40543,40603,40604,40610,40611,40612,40613,
runlist = [40667,40668,40669,40670,40671,40672,40673,40674,40675,40676,40677,40678,40679,40680,40681,40682,40683,40684,40685,40686,40687,40688,40689,40690,40691,40692,40693];%40763,40764,40765,40766];%40769,40770,40771,40772];
%pixlist = [1:61,63:70,72:85,87:100,102:112,114:123];
pixlist = 1:148;
ringlist = [1:13];
% 
% P = PLOTC('Xdata',MRA.RunData.qU,'Ydata',MRA.RunData.TBDIS,...
%     'ModelObj',MRA.ModelObj,'CovMat',MRA.FitCM,'saveplot','ON',...
%     'titleFlag','Stack','RunList',MRA.RunList,'startqU',9,...
%     'FitResult',FitResult,'RingList',ringlist);


for ri = 1:length(ringlist)
    %RA = RunAnalysis('AnaFlag','StackPixel','ring',ringlist(ri),'chi2','chi2Stat','RunNr',40667,...
    %    'fitter','minuit','exclDataStart',9,'fixPar','','ringCutFlag','');
    RA = MultiRunAnalysis('DataEffCorr','OFF','AnaFlag','SinglePixel','pixlist',pixlist(ri),'ring',ringlist(ri),'chi2','chi2Stat','RunList',runlist,...
        'fitter','minuit','exclDataStart',9,'fixPar','','pulls',[2^2,Inf,Inf,Inf,Inf,Inf]);
   % RA.RhoDScan();
   % saverhoDmin(ri) = RA.CDmin;
   % saverhoDUnc_low(ri) = RA.rhoDlowerUnc;
   % saverhoDUnc_up(ri) = RA.rhoDupperUnc;
   
   RA.FitAllSinglePixels();
   RA.PlotAllSinglePixels(pixlist);
    

  %  plotname{ri} = sprintf('plots/RhoDScan_ring%.0d.pdf',ringlist(ri));
     %RA.Fit();
     %RA.PlotFit();
     %plotname{ri} = sprintf('plots/DataSpectrumResiduals_Ring%.0d.pdf',ringlist(ri));
%     
  %export_fig(plotname{ri},'-pdf');
    close all
%     
    mnu(ri) = RA.FitResult.par(1);
    mnu_err(ri) = RA.FitResult.err(1);
    e0(ri) = RA.FitResult.par(2);
    e0_err(ri) = RA.FitResult.err(2);
    bck(ri) = RA.FitResult.par(3)*1e3;
    bck_err(ri) = RA.FitResult.err(3)*1e3;
    norm(ri) = RA.FitResult.par(4)+1;
    norm_err(ri) = RA.FitResult.err(4);
    chi2(ri) = RA.FitResult.chi2min;
    dof(ri) = RA.FitResult.dof;
% 
 RingResultsTable(ri,:) = [ringlist(ri) mnu(ri) mnu_err(ri) e0(ri) e0_err(ri) bck(ri) bck_err(ri) norm(ri) norm_err(ri) chi2(ri) dof(ri)];
    
end

% standardAveCD = 4.45e17;
% fig = figure;
% set(fig, 'Units', 'normalized', 'Position', [0.9, 0.9, 1.4, 1.5]);
% left_color = [0 0 0];
% right_color = [0 0 0];
% set(fig,'defaultAxesColorOrder',[left_color; right_color]);
% 
% yyaxis left
% errorbar(ringlist,saverhoDmin,saverhoDmin-saverhoDUnc_low,saverhoDUnc_up-saverhoDmin,...
%     's-', 'Color',rgb('CadetBlue'),'LineWidth',1.5);
% xlim([0.8 13.2])
% ylim([0 2.25*standardAveCD])
% xlabel('Ring')
% ylabel('\rho d_{min} (mol/cm^2)')
% 
% yyaxis right
% plot(ringlist,saverhoDmin./standardAveCD*100,...
%     's-', 'Color',rgb('CadetBlue'),'LineWidth',1.5);
% xlim([0.8 13.2])
% ylim([0 225])
% ylabel('100 % Col. Density = 4.45e17 mol/cm^2')
% PrettyFigureFormat;
% grid on;
% 
% plotname{ri+1} = sprintf('plots/RhoDScan_ring%.0d.pdf',ringlist(ri));
% export_fig(plotname{ri+1},'-pdf');
% 
% 
% append_pdfs('plots/RhoDScanSummaryRingNoPullsM.pdf',plotname{:})



%save('RingResultsTable.mat','RingResultsTable');
%%

titlelist = {['Neutrino mass squared [eV^2]'],['Endpoint bias [eV]'],['Background [mcps]'],['Normalization'],['Chi2 (DoF = ',num2str(RingResultsTable(1,end)),')']};

for fpd = 1:5
    figfpd = figure(fpd);
    set(figfpd,'color','white');
    FPDViewer(RingResultsTable(:,fpd*2),'ReDrawSkeleton','ON'); 
    titlehandle = title(titlelist{fpd},'interpreter','none');
    pos = get(figfpd,'position');
    set(figfpd,'position',[pos(1:2)/3.5 pos(3:4)*1.5])
    FPDName{fpd} = sprintf('plots/FPD%s.pdf',titlelist{fpd});
    export_fig([FPDName{fpd}],'-pdf');
end
close all;
    
for fpd = 1:5
    figfpd = figure(fpd);
    hold on
    bar(ringlist,RingResultsTable(:,fpd*2),'facecolor',rgb('CadetBlue'));
    if fpd < 5
        errorb(ringlist,RingResultsTable(:,fpd*2),RingResultsTable(:,fpd*2+1),'s','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black','LineWidth',1)
    end
    hold off
    titlehandle = title(titlelist{fpd},'interpreter','none');
    titlehandle.FontSize = 7;
    PrettyFigureFormat;
    BarName{fpd} = sprintf('plots/Bar%s.pdf',titlelist{fpd});
    export_fig([BarName{fpd}],'-pdf');
end
close all;
delete('plots/ResultsSummaryRing.pdf');
append_pdfs('plots/ResultsSummaryRingJustDT.pdf',FPDName{:},BarName{:},plotname{:})
    
% pixlist40604 = [1:100,102,104:112,115,116,118:123];

