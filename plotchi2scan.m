function plotchi2scan(savename)
    load(savename);
    
    if contains(savename,'2D')
        [~,s]=contour(etaScanPoints,mnuScanPoints,Chi2,50);
        %s.FaceColor='interp';
        s.LineWidth = 2;
        xlabel('\eta','FontSize',12);
        ylabel('m_{\nu}^{2} (eV^{2})','FontSize',12);
        zlabel('\chi^2','FontSize',12);
        PRLFormat;
        set(gca, 'YScale', 'log');
        set(gca, 'XScale', 'log');
        colorbar;
        hold on;
        [~,t]=contour(etaScanPoints,mnuScanPoints,Chi2,[0 4.61]);
        t.LineWidth = 2;
        t.ShowText = 'ON';
        PrettyFigureFormat;
        hold off;
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