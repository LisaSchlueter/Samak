% Understanding systematics breakdown
    DataType = 'Twin';
    SavePlt = 'ON';
    savedir = [getenv('SamakPath'),'ksn2ana/ksn2_Systematics/results/'];
    savename = sprintf('%sksn2_SystBreakdown_StatOverSyst_%s_%.0feV_RasterScan%s.mat',savedir,DataType,40,'ON');
    if exist(savename,'file')
        d = importdata(savename);
        fprintf('load file %s \n',savename);
    else
        return    
    end
    plotdir = strrep(savedir,'results','plots');
    %%  
    f1 = figure('Units','normalized','Position',[-0.1,0.1,0.5,0.6]);
    s1 = subplot(4,1,1);
    x = -sqrt(d.mNu4Sq(1,:));
    Idx = 1; %'LongPlasma','BkgPT' (2,3)
    pStat = plot(x,d.sin2t4_Stat(1,:).^2./d.sin2t4_Tot(1,:).^2,'LineWidth',3,'Color',rgb('Silver'),'LineStyle','-');
    hold on;
    pTot = plot(x,d.sin2t4_Sys(Idx,:).^2./d.sin2t4_Tot(1,:).^2,'-k','LineWidth',3);
    xlabel(sprintf('- {\\itm}_4^ (eV)'));
    ylabel(sprintf('\\sigma^2 / \\sigma_{total}^2'));
    leg = legend('Stat. only','All syst. effects','Location','west');
    PrettyLegendFormat(leg);
    ylim([0 1])
    ax1 = gca;
    PrettyFigureFormat('FontSize',20);
    
    s2 = subplot(4,1,2);
    bar(d.qU-18574,d.TimeSubRun./(60*60),'FaceColor',rgb('DodgerBlue'),'EdgeColor','none');
   % xlabel(sprintf('{\\itqU} - 18574 (eV)'))
    ylabel(sprintf('Time (h)'));
    ylim([0 70]);
    PrettyLegendFormat(leg);
    ax2 = gca;
    xticklabels('');
    PrettyFigureFormat('FontSize',20);
    
    s3 = subplot(4,1,3);
    plot(d.qU-18574,ones(numel(d.qU),1),'k-','LineWidth',1.5);
    hold on;
    plot(d.qU-18574,d.TBDIS_Signal./d.BkgCounts,'-','MarkerSize',18,'LineWidth',2.5);
    %xlabel(sprintf('{\\itqU} - 18574 (eV)'))
    ylabel(sprintf('{\\itS} / {\\itB}'));
    set(gca,'YScale','log');
    ylim([0.005 500])
    yticks([0.01 1 100])
    xticklabels('');
    ax3 = gca;
    PrettyFigureFormat('FontSize',20);
      
    s4 = subplot(4,1,4);
    plot(d.qU-18574,1e-05.*(d.TBDIS_Signal+d.BkgCounts),'.-','LineWidth',2.5,'Color',rgb('DodgerBlue'),'MarkerSize',18);
    hold on;
    plot(d.qU-18574,d.BkgCounts.*1e-05,':','LineWidth',2,'Color',rgb('ForestGreen'));
    xlabel(sprintf('{\\itqU} - 18574 (eV)'))
    ylabel(sprintf('Counts (10^5)'))
    ylim([-0.5 5.7]);
    yticks([1 3 5]);
    leg = legend(sprintf('{\\itS} + {\\itB}'),sprintf('{\\itB}'),'Location','northeast');
    leg.NumColumns = 2;
    PrettyLegendFormat(leg);
    ax4 = gca;
     PrettyFigureFormat('FontSize',20);
    
    linkaxes([s1,s2,s3,s4],'x');
    xlim([-37,0]);
    
    ax1.Position(4) = 0.16;
    ax2.Position(4) = ax1.Position(4);
    ax3.Position(4) = ax1.Position(4);
    ax4.Position(4) = ax1.Position(4);
    ax1.Position(2) = 0.82;
    ax2.Position(2) = 0.54;
    ax3.Position(2) = 0.34;
    ax4.Position(2) = 0.14;
    ax1.YLabel.Position(1) = -39.5;
    ax2.YLabel.Position(1) = ax1.YLabel.Position(1);
    ax3.YLabel.Position(1) = ax1.YLabel.Position(1);
    ax4.YLabel.Position(1) = ax1.YLabel.Position(1);
    %
    if strcmp(SavePlt,'ON')
        plotname = sprintf('%sksn2_SysBreakdown_MTD_%s_RasterScan%s.png',plotdir,DataType,'ON');
        print(f1,plotname,'-dpng','-r400');
        fprintf('save plot to %s \n',plotname);
       % export_fig(f1,strrep(plotname,'png','pdf'));
    end
%%       alternative plot

f2 = figure('Units','normalized','Position',[-0.1,0.1,0.6,0.6]);
s1 = subplot(2,3,[1,2,4,5]);
pstat =plot(d.sin2t4_Stat(1,:),d.mNu4Sq(1,:),'Color',rgb('Silver'),'LineWidth',3);
hold on;
ptot =plot(d.sin2t4_Tot(1,:),d.mNu4Sq(1,:),'-.','Color',rgb('Black'),'LineWidth',3);
leg = legend([pstat,ptot],'Stat. only','Total (stat. and syst.)','Location','northeast'); 
PrettyLegendFormat(leg); legend boxoff
set(gca,'XScale','log')
set(gca,'YScale','log')
PrettyFigureFormat('FontSize',22);
ax1 = gca;
xlabel(sprintf('|{\\itU}_{e4}|^2'));
ylabel(sprintf('{\\itm}_4^2 (eV)'));
xlim([5e-03 0.5]);

s2 = subplot(2,3,[3,6]);
Colors = colormap('jet');
pTot = plot(d.sin2t4_Sys(1,:).^2./d.sin2t4_Tot(1,:).^2,d.mNu4Sq(1,:),'LineWidth',3,'Color',rgb('Black'));
hold on;
pStat= plot(d.sin2t4_Stat(1,:).^2./d.sin2t4_Tot(1,:).^2,d.mNu4Sq(1,:),'LineWidth',3,'Color',rgb('Silver'),'LineStyle','-');
% pPlasma = plot(d.sin2t4_Sys(2,:).^2./d.sin2t4_Tot(1,:).^2,d.mNu4Sq(2,:),'LineWidth',1.5,'Color',Colors(floor((2)*256/(d.nSys)),:),'LineStyle','-.');
% pPenning = plot(d.sin2t4_Sys(3,:).^2./d.sin2t4_Tot(1,:).^2,d.mNu4Sq(3,:),'LineWidth',1.5,'Color',Colors(floor((3)*256/(d.nSys)),:),'LineStyle',':');
% pNonPois = plot(d.sin2t4_Sys(4,:).^2./d.sin2t4_Tot(1,:).^2,d.mNu4Sq(4,:),'LineWidth',1.5,'Color',Colors(floor((4)*256/(d.nSys)),:),'LineStyle','--');
%leg = legend([pStat,pTot,pPlasma,pPenning,pNonPois],'Stat. only','All syst.','Plasma','Penning bkg','Non-Pois. bkg','Location','northeast');
leg = legend([pStat,pTot],'Stat. only','All syst.','Location','north');

PrettyLegendFormat(leg);
PrettyFigureFormat('FontSize',22);
set(gca,'YScale','log');
ax2 = gca;
ax2.YAxisLocation = 'right';
xlabel(sprintf('\\sigma^2 / \\sigma_{total}^2'));
xlim([0 1]);
linkaxes([s1,s2],'y');
ax1.Position(3) = 0.55;
ax1.Position(2) = 0.15;
ax2.Position(2) = ax1.Position(2);
ax1.Position(1) = 0.1;
ax2.Position(1) = 0.66;


if strcmp(SavePlt,'ON')
    plotname = sprintf('%sksn2_SysBreakdown_RatioContour_%s_RasterScan%s.png',plotdir,DataType,'ON');
    print(f2,plotname,'-dpng','-r400');
    fprintf('save plot to %s \n',plotname);
    %export_fig(f2,strrep(plotname,'png','pdf'));
end