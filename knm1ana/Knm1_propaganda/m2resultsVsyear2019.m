function m2resultsVsyear2019()
% Plot of Squared neutrino mass values obtained from tritium beta decay
% in the  period 1990-2019 plotted against the year of publication 
% Previous: results from the more recent experiments in Mainz and Troitsk
% New: results from KATRIN KNM1

%% Load Data
[mBetaSquared , mBetaSquaredPDG]  = mBetaSquaredExperiments();

%% Plot
maintitle=sprintf('Squared neutrino mass values obtained from tritium \\beta -decay in the period 1990-2019');
savefile=sprintf('plots/m2resultsVsyear2019.png');
fig1 = figure('Name','Test','NumberTitle','off','rend','painters','pos',[10 10 1400 600]);
a=annotation('textbox', [0 0.9 1 0.1], ...
    'String', maintitle, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
a.FontSize=26;a.FontWeight='normal';
for i=1:12
    h(i)=errorbar(mBetaSquared(i).Year, mBetaSquared(i).mBetaSquared,mBetaSquared(i).Tot,...
    'Marker',mBetaSquared(i).Marker,'MarkerSize',12,'LineWidth',1,'Color','Black',...
    'MarkerEdgeColor',mBetaSquared(i).Color,'MarkerFaceColor',mBetaSquared(i).Color);
    h(i).CapSize = 0; h(i).LineStyle= 'none';h(i).LineWidth= 3;
    hold on
end
[pdg1 pdg2]= boundedline([mBetaSquared(1).Year-1 mBetaSquared(end).Year+1],[-0.6 -0.6],1.9,'alpha','cmap',rgb('CadetBlue'));
% Zoom KATRIN
% zoomPlot to highlight a portion of the major plot
l=legend(h,[mBetaSquared(1).Experiment ' ' mBetaSquared(2).Reference],...
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
    [mBetaSquared(12).Experiment ' ' mBetaSquared(12).Reference],...
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
xlim([mBetaSquared(1).Year-1 mBetaSquared(12).Year+1]);
% Zoom
hold on
%rectangle('Position',[2000 -225 15 175]);
z = axes(gcf,'position',[0.37 0.17 .52 .41],'FontSize',16,'XAxisLocation', 'top');
box on 
for i=10:13
    hz(i)=errorbar(mBetaSquared(i).Year, mBetaSquared(i).mBetaSquared,mBetaSquared(i).Tot,...
    'Marker',mBetaSquared(i).Marker,'MarkerSize',12,'LineWidth',1,'Color','Black',...
    'MarkerEdgeColor',mBetaSquared(i).Color,'MarkerFaceColor',mBetaSquared(i).Color);
    hz(i).MarkerSize = 16; hz(i).CapSize = 0; z.FontSize= 20;
    hz(i).LineStyle= 'none';hz(i).LineWidth= 5;legend hide
    xlim([2004 2021])
hold on
end
xticks([2005 2011 mBetaSquared(12).Year-0.5]);
yticks([-3 -2 -1 0 1 2 ]);
[pdg1 pdg2]=boundedline([mBetaSquared(1).Year-1 mBetaSquared(end).Year+1],[-0.6 -0.6],1.9,'alpha','cmap',rgb('CadetBlue'));
set(gca,'XAxisLocation', 'top')
grid on
ylabel('m^2(\nu_e) c^4 (eV^2)');
l=legend([ hz(10) hz(11)  hz(12) pdg2],...
    [mBetaSquared(10).Experiment ' '],...
    [mBetaSquared(11).Experiment ' '],...
    [mBetaSquared(12).Experiment ' '],...
    [mBetaSquaredPDG.Label ' '],...
    'FontSize',20,'Location','SouthEast');
l.NumColumns=2;
legend boxoff;
export_fig(gcf,savefile,'-r300');