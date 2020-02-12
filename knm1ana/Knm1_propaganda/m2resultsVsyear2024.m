function m2resultsVsyear2024()
% Plot of Squared neutrino mass values obtained from tritium beta decay
% in the  period 1990-2019 plotted against the year of publication 
% Previous: results from the more recent experiments in Mainz and Troitsk
% New: results from KATRIN KNM1

%% Load Data
[mBetaSquared , mBetaSquaredPDG] = mBetaSquaredExperiments();

%% Plot
maintitle=sprintf('Squared neutrino mass values obtained from tritium \\beta -decay in the period 1990-2019');
savefile=sprintf('plots/m2resultsVsyear2024.png');
fig1 = figure('Name','Test','NumberTitle','off','rend','painters','pos',[10 10 1400 600]);
a=annotation('textbox', [0 0.9 1 0.1], ...
    'String', maintitle, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
a.FontSize=26;a.FontWeight='normal';
for i=1:13
    h(i)=errorbar(mBetaSquared(i).Year, mBetaSquared(i).mBetaSquared,mBetaSquared(i).Tot,...
    'Marker',mBetaSquared(i).Marker,'MarkerSize',12,'LineWidth',1,'Color','Black',...
    'MarkerEdgeColor',mBetaSquared(i).Color,'MarkerFaceColor',mBetaSquared(i).Color);
    h(i).CapSize = 0; h(i).LineStyle= 'none';h(i).LineWidth= 3;
    hold on
end
line([mBetaSquared(1).Year-1 mBetaSquared(13).Year+1],[0 0 ],'LineStyle','--','LineWidth',1,'Color',rgb('SteelBlue'));
% Zoom KATRIN
% zoomPlot to highlight a portion of the major plot
l=legend(h,[mBetaSquared(1).Experiment ' ' mBetaSquared(1).Reference],...
    [mBetaSquared(2).Experiment ' ' mBetaSquared(2).Reference],...
       [mBetaSquared(3).Experiment ' ' mBetaSquared(3).Reference],...
    [mBetaSquared(4).Experiment ' ' mBetaSquared(4).Reference],...
    [mBetaSquared(5).Experiment ' ' mBetaSquared(5).Reference],...
    [mBetaSquared(6).Experiment ' ' mBetaSquared(6).Reference],...
    [mBetaSquared(7).Experiment ' ' mBetaSquared(7).Reference],...
    [mBetaSquared(8).Experiment ' ' mBetaSquared(8).Reference],...
    [mBetaSquared(9).Experiment ' ' mBetaSquared(9).Reference],...
    [mBetaSquared(10).Experiment ' ' mBetaSquared(10).Reference],...
    [mBetaSquared(11).Experiment ' ' mBetaSquared(11).Reference],...
    [mBetaSquared(13).Experiment ' ' mBetaSquared(13).Reference],...
    'FontSize',12,'Location','north');
l.NumColumns=4;
legend('boxoff');
xlabel('year');
ylabel('m^2(\nu_e) c^4 (eV^2)');
grid on
FontName = 'Arial';
FontSize = 22;
set(gca,'FontName',FontName,'FontSize',FontSize);
set(get(gca,'XLabel'),'FontName',FontName,'FontSize',FontSize);
set(get(gca,'YLabel'),'FontName',FontName,'FontSize',FontSize);
set(gca,'TickLength',[.01 .01]);
set(gca,'XMinorTick','off');
set(gca,'YMinorTick','on');
set(gcf,'Color','w');
box on;
xlim([mBetaSquared(1).Year-1 mBetaSquared(13).Year+1]);
% Zoom
hold on
z = axes(gcf,'position',[0.37 0.2 .5 .38],'FontSize',16,'XAxisLocation', 'top');
box on 
for i=12:13
    hz(i)=errorbar(mBetaSquared(i).Year, mBetaSquared(i).mBetaSquared,mBetaSquared(i).Tot,...
    'Marker',mBetaSquared(i).Marker,'MarkerSize',12,'LineWidth',1,'Color','Black',...
    'MarkerEdgeColor',mBetaSquared(i).Color,'MarkerFaceColor',mBetaSquared(i).Color);
    hz(i).MarkerSize = 14; hz(i).CapSize = 0; z.FontSize= 20;
    hz(i).LineStyle= 'none';hz(i).LineWidth= 5;legend hide
    xlim([2019 2025])
hold on
end
xticks([2019.5 mBetaSquared(13).Year]);
yticks([-1 1]);
ylim([-2.1 1.1]);
line([mBetaSquared(1).Year-1 mBetaSquared(13).Year+1],[0 0 ],'LineStyle','--','LineWidth',1,'Color',rgb('SteelBlue'));
set(gca,'XAxisLocation', 'top')
grid on
ylabel('m_\beta c^4 (eV^2)');
l=legend([hz(12) hz(13)],...
    [mBetaSquared(12).Label sprintf(' - 1 \\sigma = %.1f',mBetaSquared(12).Stat) ' (stat) + ' sprintf('%.1f',mBetaSquared(12).Sys)  ' (sys) eV^2'],...
    [mBetaSquared(13).Label sprintf(' - 1 \\sigma = %.3f',mBetaSquared(13).Stat) ' (stat) + ' sprintf('%.3f',mBetaSquared(13).Sys)  ' (sys) eV^2'],...
    'FontSize',20,'Location','South');
legend boxoff;
export_fig(gcf,savefile,'-r300');