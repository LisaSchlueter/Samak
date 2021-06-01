function plotchi2scan(savename)
    
    if exist(savename,'file')
        load(savename);
    end
    
    if contains(savename,'2D')
        fig1=figure(1);
        contour(mnuScanPoints,etaScanPoints,(Chi2-GlobalChi2Min).',[0 4.61]);
         c=colorbar;
         c.Label.String = '\Delta\chi^2';
%         title('\chi^2');
         %s.LineWidth = 2;
         PRLFormat;
         %set(gca, 'YScale', 'log');
         %set(gca, 'XScale', 'log');
        hold on;
         s=imagesc(mnuScanPoints,etaScanPoints,(Chi2-GlobalChi2Min).',[0 8]);
         %s.FaceColor='interp';
         [~,v]=contour(mnuScanPoints,etaScanPoints,(Chi2-GlobalChi2Min).',[0 4.61]);
        %[~,t]=contour(mnuScanPoints,etaScanPoints,(Chi2-GlobalChi2Min).',[0 2.3]);
        %t.LineWidth = 2;
        %t.LineColor = 'White';
        %t.ShowText = 'ON';
        %[~,u]=contour(mnuScanPoints,etaScanPoints,(Chi2-GlobalChi2Min).',[0 6.18]);
        %u.LineWidth = 2;
        %u.LineColor = 'White';
        %u.ShowText = 'ON';
        v.LineWidth = 2;
        v.LineColor = 'White';
        v.ShowText = 'ON';
        ylabel('\eta','FontSize',12);
        xlabel('m_{\nu}^{2} (eV^{2})','FontSize',12);
        zlabel('\Delta\chi^2','FontSize',12);
        %title('\Delta\chi^2');
        %text(0,2.5e11,'1\sigma C.L.','FontSize',14);
        %text(-0.5,4e11,'2\sigma C.L.','FontSize',14);
        %text(-0.25,6e11,'3\sigma C.L.','FontSize',14);
        xlim([0 2.3]);
        ylim([0 inf]);
        PrettyFigureFormat;
        hold off;
        
        fig2=figure(2);
        imagesc(etaScanPoints,mnuScanPoints,E0);
        colorbar;
        title('E_0');
        xlabel('\eta','FontSize',12);
        ylabel('m_{\nu}^{2} (eV^{2})','FontSize',12);
        PRLFormat;
        
        fig3=figure(3);
        imagesc(etaScanPoints,mnuScanPoints,Norm);
        colorbar;
        title('Normalization');
        xlabel('\eta','FontSize',12);
        ylabel('m_{\nu}^{2} (eV^{2})','FontSize',12);
        PRLFormat;
        
        fig4=figure(4);
        imagesc(etaScanPoints,mnuScanPoints,Bkg);
        colorbar;
        title('Background');
        xlabel('\eta','FontSize',12);
        ylabel('m_{\nu}^{2} (eV^{2})','FontSize',12);
        PRLFormat;
        
    elseif contains(savename,'EtaFit')
        fig1=figure(10);
        boundedline(sqrt(mNu),etaFit,[(eta-etaFit); eta-etaFit]');
        xlabel('m_{\nu} (eV)');
        ylabel('\eta');
        legend('\eta best fit with 90% C.L.','box','off','Location','northwest');
        PRLFormat;
        
        FitStyleArg = {'o','Color','k','LineWidth',1.0,'MarkerFaceColor',rgb('Black'),'MarkerSize',4,'MarkerEdgeColor',rgb('Black')};

        fig2=figure(20);
        pE0=errorbar(mNu,E0,E0_err,FitStyleArg{:},'CapSize',0);
        hold on;
        fitobject = fit(permute(mNu,[2 1]),permute(E0,[2 1]),'poly1');
        plot(fitobject,mNu,E0);
        xlabel('m_{\nu}^2 (ev^2)','FontSize',12);
        ylabel('E_0^{fit} (eV)','FontSize',12);
        lgd=legend([pE0],{'E_0^{fit} with linear trend'},'Location','northwest','box','off');
        lgd.FontSize = 12;
        PrettyFigureFormat;
        hold off;
        
        fig3=figure(30);
        pChi2=plot(mNu,Chi2,'LineWidth',2);
        xlabel('m_{\nu}^2 (ev^2)','FontSize',12);
        ylabel('\chi^2','FontSize',12);
        lgd=legend([pChi2],{'\chi^2'},'Location','northwest','box','off');
        lgd.FontSize = 12;
        PrettyFigureFormat;
        hold off;
        
        fig4=figure(40);
        pNorm=errorbar(mNu,Norm,Norm_err,FitStyleArg{:},'CapSize',0);
        hold on;
        fitobject = fit(permute(mNu,[2 1]),permute(Norm,[2 1]),'poly1');
        plot(fitobject,mNu,Norm);
        xlabel('m_{\nu}^2 (ev^2)','FontSize',12);
        ylabel('N','FontSize',12);
        lgd=legend([pNorm],{'Normalization with linear trend'},'Location','northwest','box','off');
        lgd.FontSize = 12;
        PrettyFigureFormat;
        hold off;
        
        fig5=figure(50);
        pBkg=errorbar(mNu,Bkg,Bkg_err,FitStyleArg{:},'CapSize',0);
        hold on;
        fitobject = fit(permute(mNu,[2 1]),permute(Bkg,[2 1]),'poly1');
        plot(fitobject,mNu,Bkg);
        xlabel('m_{\nu}^2 (ev^2)','FontSize',12);
        ylabel('Bkg (mcps)','FontSize',12);
        lgd=legend([pBkg],{'Background rate with linear trend'},'Location','northwest','box','off');
        lgd.FontSize = 12;
        PrettyFigureFormat;
        hold off;
        
    elseif strcmp(savename,'KRN1Raster')
        load('./RelicNuBkg/UpperLimits/EtaFit_mNu0_1_SystON_DataTypeRealSinglemnuSq_TwinBias_mnuSq0_E0 Norm Bkg etaDeltaChi2_1.mat');
        eta90=eta;
        etaFit90=etaFit;
        load('./RelicNuBkg/UpperLimits/EtaFit_mNu0_1_SystON_DataTypeRealSinglemnuSq_TwinBias_mnuSq0_E0 Norm Bkg etaDeltaChi2_4.mat');
        eta95=eta;
        etaFit95=etaFit;
        load('./RelicNuBkg/UpperLimits/EtaFit_mNu0_1_SystON_DataTypeRealSinglemnuSq_TwinBias_mnuSq0_E0 Norm Bkg etaDeltaChi2_9.mat');
        hold on;
        [~,b] = boundedline(sqrt(mNu),etaFit90,eta90,'cmap',[0 0.75 0.75]);
        [~,d] = boundedline(sqrt(mNu),etaFit95,eta95,'cmap',[0 0.5 0.5],'alpha');
        [~,f] = boundedline(sqrt(mNu),etaFit,eta,'cmap',[0 0.25 0.25],'alpha');
        %legend([b d f],'\eta 90% C.L.','\eta 95% C.L.','\eta 99% C.L.');
        text(0.1,3.3e11,'1\sigma C.L.','FontSize',14);
        text(0.4,5e11,'2\sigma C.L.','FontSize',14);
        text(0.7,6.5e11,'3\sigma C.L.','FontSize',14);
        ylim([0 inf]);
        xlabel('m_{\nu}');
        ylabel('\eta');
        PrettyFigureFormat;
    else

        eta=(0:(Netabins-1))*((etafactor*10^(etarange))/(Netabins-1));

        fig1=figure(1);
        plot(eta,Chi2,'LineWidth',2);
        xlabel('\eta','FontSize',12);
        ylabel('\chi^2','FontSize',12);
        PrettyFigureFormat;

        FitStyleArg = {'o','Color','k','LineWidth',1.0,'MarkerFaceColor',rgb('Black'),'MarkerSize',4,'MarkerEdgeColor',rgb('Black')};

        fig2=figure(2);
        pmnuSq=errorbar(eta,mnuSq,mnuSq_err,FitStyleArg{:},'CapSize',0);
        hold on;
        fitobject = fit(permute(eta,[2 1]),permute(mnuSq,[2 1]),'poly1');
        plot(fitobject,eta,mnuSq);
        xlabel('\eta','FontSize',12);
        ylabel('m_{\nu}^{2} (eV^{2})','FontSize',12);
        xlim([eta(1)-0.1*eta(end) eta(end)+0.1*eta(end)]);
        lgd=legend([pmnuSq],{'m_{\nu}^{2} with linear trend'},'Location','northwest','box','off');
        lgd.FontSize = 12;
        PrettyFigureFormat;
        hold off;

        fig3=figure(3);
        pE0=errorbar(eta,E0,E0_err,FitStyleArg{:},'CapSize',0);
        hold on;
        fitobject = fit(permute(eta,[2 1]),permute(E0,[2 1]),'poly1');
        plot(fitobject,eta,E0);
        xlabel('\eta','FontSize',12);
        ylabel('E_{0}^{fit} (eV)','FontSize',12);
        xlim([eta(1)-0.1*eta(end) eta(end)+0.1*eta(end)]);
        lgd=legend([pE0],{'E_{0}^{fit} with linear trend'},'Location','northwest','box','off');
        lgd.FontSize = 12;
        PrettyFigureFormat;
        hold off;

        fig4=figure(4);
        pB=errorbar(eta,Bkg,Bkg_err,FitStyleArg{:},'CapSize',0);
        hold on;
        fitobject = fit(permute(eta,[2 1]),permute(Bkg,[2 1]),'poly1');
        plot(fitobject,eta,Bkg);
        xlabel('\eta','FontSize',12);
        ylabel('Background (cps)','FontSize',12);
        xlim([eta(1)-0.1*eta(end) eta(end)+0.1*eta(end)]);
        lgd=legend([pB],{'Bkg rate with linear trend'},'Location','northeast','box','off');
        lgd.FontSize = 12;
        PrettyFigureFormat;
        hold off;

        fig5=figure(5);
        pN=errorbar(eta,Norm,Norm_err,FitStyleArg{:},'CapSize',0);
        hold on;
        fitobject = fit(permute(eta,[2 1]),permute(Norm,[2 1]),'poly1');
        plot(fitobject,eta,Norm);
        xlabel('\eta','FontSize',12);
        ylabel('N','FontSize',12);
        xlim([eta(1)-0.1*eta(end) eta(end)+0.1*eta(end)]);
        lgd=legend([pN],{'Normalization with linear trend'},'Location','northeast','box','off');
        lgd.FontSize = 12;
        PrettyFigureFormat;
        hold off;     
    end
end