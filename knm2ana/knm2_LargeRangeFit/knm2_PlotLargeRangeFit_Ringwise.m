[~, E0, E0Err,FitResults,Rate_m,Rate_d,qU_d,TimeSecSubRun] = knm2_LargeRangeFit_Ringwise;

Colors = {rgb('DodgerBlue');rgb('Orange');rgb('FireBrick');rgb('ForestGreen')};
f81 = figure('units','normalized','pos',[0.1, 0.1,0.5,0.8]);
for i=1:3
    qUmean = squeeze(mean(qU_d(i,:,:),2));
    s1 = subplot(3,1,i);
    %pFit = plot(qUmean-18574,squeeze(Rate_m(i,:,:))','-','LineWidth',1.5,...
    %    'Color',Colors);
    hold on;
    for r=1:4
        pFit = plot(qUmean-18574,squeeze(Rate_m(i,r,:))','-','LineWidth',1.5,...
            'Color',[Colors{r}]);
        hold on;
        e1= errorbar(qUmean-18574,squeeze(Rate_d(i,r,:))',sqrt(squeeze(Rate_d(i,r,:)))','.',...
            'CapSize',0,'MarkerSize',20,'Color',pFit.Color);
    end
    xlim([-91 -43]);
    hold off;
    PrettyFigureFormat('FontSize',20)
    leg = legend([e1,pFit(end)],sprintf('KATRIN data with 1\\sigma error bars'),...
        'Fit result','EdgeColor',rgb('Silver'),'FontSize',get(gca,'FontSize'));
    leg.Title.String = sprintf('Period %.0f',i);
    leg.Title.FontWeight = 'normal';
    leg.LineWidth = 1;
  
    if i==1
        title(sprintf('4 Pseudo-Rings , fit range [-90 to -45, +135] eV'),'FontWeight','normal');
    elseif i==2
        ax = gca;
        mypos = ax.Position;
        ax.Position = [mypos(1) mypos(2)+0.04 mypos(3:4)];
    elseif i==3
         ax2 = gca;
         mypos = ax2.Position;
         ax2.Position = [mypos(1) mypos(2)+0.08 mypos(3:4)];
         xlabel('Retarding potential - 18574 (eV)');
    end
    ylabel(sprintf('Rate (cps)'))
    ylim([-20,1.15*max(max(Rate_m(1,:,:)))])
end

plotdir = [getenv('SamakPath'),'knm2ana/knm2_LargeRangeFit/plots/'];
MakeDir(plotdir);
plotname = sprintf('%sknm2_LargeRangeFit_Ringwise.pdf',plotdir);
export_fig(plotname);
fprintf('save plot to %s \n',plotname);


%% plot residuals
f82 = figure('units','normalized','pos',[0.1, 0.1,0.5,0.8]);
for i=1:3
    qUmean = squeeze(mean(qU_d(i,:,:),2));
       Residuals = squeeze((Rate_d(i,:,:)-Rate_m(i,:,:))./sqrt(Rate_d(i,:,:)).*sqrt(TimeSecSubRun(i,:,:)));
       s1 = subplot(3,1,i);
       pref = plot(linspace(-100,-30,10),zeros(10,1),'-','LineWidth',1.5,'Color',rgb('DimGray'));
       hold on;
       for r=1:4
           pFit = plot(qUmean-18574,Residuals(r,:),'.-.','LineWidth',1.5,'MarkerSize',20,'Color',Colors{r});
       end
    xlim([-91 -43]);
    hold off;
    PrettyFigureFormat('FontSize',20)
    if i==1
        title(sprintf('4 Pseudo-Rings , fit range [-90 to -45, +135] eV'),'FontWeight','normal');
    elseif i==2
        ax = gca;
        mypos = ax.Position;
        ax.Position = [mypos(1) mypos(2)+0.04 mypos(3:4)];
    elseif i==3
         ax2 = gca;
         mypos = ax2.Position;
         ax2.Position = [mypos(1) mypos(2)+0.08 mypos(3:4)];
         xlabel('Retarding potential - 18574 (eV)');
    end
    ylabel(sprintf('Residuals (\\sigma)'))
    ylim([-3 3])
    
end
plotnameRes = strrep(plotname,'.pdf','_Residuals.pdf');
export_fig(plotnameRes);
fprintf('save plot to %s \n',plotnameRes);

