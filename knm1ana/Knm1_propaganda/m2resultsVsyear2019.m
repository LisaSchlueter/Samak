function [leg,ax1,z,b1] = m2resultsVsyear2019(varargin)
p=inputParser;
p.addParameter('Combi','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('LegPos','north',@(x)ischar(x));
p.addParameter('LegOpt','Full',@(x)ismember(x,{'Full','Short'}));
p.addParameter('SavePlot','OFF',@(x)ismember(x,{'ON','OFF'}));

p.parse(varargin{:});

Combi    = p.Results.Combi;
LegPos   = p.Results.LegPos;
LegOpt   = p.Results.LegOpt;
SavePlot = p.Results.SavePlot;

% Plot of Squared neutrino mass values obtained from tritium beta decay
% in the  period 1990-2019 plotted against the year of publication 
% Previous: results from the more recent experiments in Mainz and Troitsk
% New: results from KATRIN KNM1
if strcmp(Combi,'OFF')
    LocalFontSize = 32;
    MarkerSize = 12;
else
    LocalFontSize = 20;
    MarkerSize = 8;
end
%% Load Data
[mBetaSquared , mBetaSquaredPDG]  = mBetaSquaredExperiments();

%% Plot
maintitle=sprintf('Squared neutrino mass values obtained from tritium \\beta -decay in the period 1990-2019');
savefile=sprintf('plots/m2resultsVsyear2019.png');
if strcmp(Combi,'OFF')
fig1 = figure('Name','Test','NumberTitle','off','pos',[10 10 1400 600]);
end
% a=annotation('textbox', [0 0.9 1 0.1], ...
%     'String', maintitle, ...
%     'EdgeColor', 'none', ...
%     'HorizontalAlignment', 'center');
% a.FontSize=26;a.FontWeight='normal';
[pdg1 pdg2]= boundedline([mBetaSquared(1).Year-1 mBetaSquared(end).Year+1],[-0.6 -0.6],1.9,'alpha','cmap',rgb('CadetBlue'));
hold on
for i=1:12
    h(i)=errorbar(mBetaSquared(i).Year, mBetaSquared(i).mBetaSquared,mBetaSquared(i).Tot,...
    'Marker',mBetaSquared(i).Marker,'MarkerSize',MarkerSize,'LineWidth',1,'Color','Black',...
    'MarkerEdgeColor',mBetaSquared(i).Color,'MarkerFaceColor',mBetaSquared(i).Color);
    h(i).CapSize = 0; h(i).LineStyle= 'none';h(i).LineWidth= 2;
end
set(gca,'TickDir','out');
% remove top and right ticks
ax1 = gca;
mypos1 = ax1.Position;
set(ax1,'box','off','color','none')% set box property to off and remove background color
b1 = axes('Position',[mypos1(1) mypos1(2) mypos1(3) mypos1(4)],...
    'box','on','xtick',[],'ytick',[],'LineWidth',1.5);% create new, empty axes with box but without ticks
axes(ax1)% set original axes as active
linkaxes([ax1 b1]) % link axes in case of zooming

switch LegOpt
    case 'Full'
        leg=legend(h,[mBetaSquared(1).Experiment ' ' mBetaSquared(2).Reference],...
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
            'FontSize',13.5,'Location',LegPos);
    case 'Short'
        leg=legend(h,[mBetaSquared(1).Experiment ' ' sprintf('(%.0f)',mBetaSquared(1).Year)],...
            [mBetaSquared(2).Experiment ' ' sprintf('(%.0f)',mBetaSquared(2).Year)],...
            [mBetaSquared(3).Experiment ' ' sprintf('(%.0f)',mBetaSquared(3).Year)],...
            [mBetaSquared(4).Experiment ' ' sprintf('(%.0f)',mBetaSquared(4).Year)],...
            [mBetaSquared(5).Experiment ' ' sprintf('(%.0f)',mBetaSquared(5).Year)],...
            [mBetaSquared(6).Experiment ' ' sprintf('(%.0f)',mBetaSquared(6).Year)],...
            [mBetaSquared(7).Experiment ' ' sprintf('(%.0f)',mBetaSquared(7).Year)],...
            [mBetaSquared(8).Experiment ' ' sprintf('(%.0f)',mBetaSquared(8).Year)],...
            [mBetaSquared(9).Experiment ' ' sprintf('(%.0f)',mBetaSquared(9).Year)],...
            [mBetaSquared(10).Experiment ' ' sprintf('(%.0f)',mBetaSquared(10).Year)],...
            [mBetaSquared(11).Experiment ' ' sprintf('(%.0f)',mBetaSquared(11).Year)],...
            ['KATRIN' ' ' sprintf('(%.0f)',mBetaSquared(12).Year)],...
            'FontSize',13.5,'Location',LegPos);
end

if strcmp(Combi,'OFF')
leg.NumColumns=4;
else
    leg.NumColumns=4;
    leg.FontSize = 18;
end
legend('boxoff');
xlabel('Year');
ylabel(sprintf('{\\itm}_\\nu^2 (eV^2)'));%c^4 grid on
ax1 = gca;
xticks([]);
%FontName = 'Arial';

%set(gca,'FontName',FontName,'FontSize',FontSize);
%set(get(gca,'XLabel'),'FontName',FontName,'FontSize',LocalFontSize);
%set(get(gca,'YLabel'),'FontName',FontName,'FontSize',LocalFontSize);
%set(gca,'TickLength',[.01 .01]);
%%set(gcf,'Color','w');
%set(gca,'YMinorTick','on');

PRLFormat
box off
set(gca,'FontSize',LocalFontSize);
set(get(gca,'XLabel'),'FontSize',LocalFontSize+4);
set(get(gca,'YLabel'),'FontSize',LocalFontSize+4);
set(gca,'XMinorTick','off');
%box on;
ylim([-250 110])
xlim([mBetaSquared(1).Year-1 mBetaSquared(12).Year+1]);
% Zoom
hold on

% Zoom KATRIN
% zoomPlot to highlight a portion of the major plot
%rectangle('Position',[2000 -225 15 175]);
if strcmp(Combi,'OFF')
    z = axes(gcf,'position',[0.37 0.17 .52 .41],'FontSize',LocalFontSize,'XAxisLocation', 'top');
else
    z = axes(gcf,'position',[0.4 0.65 .4 .15],'FontSize',LocalFontSize-4,'XAxisLocation', 'top');
    
end
box on 
[pdg1 pdg2]=boundedline([mBetaSquared(1).Year-1 mBetaSquared(end).Year+1],[-0.6 -0.6],1.9,...
    'alpha','cmap',rgb('CadetBlue'));
hold on;
for i=10:13
    hz(i)=errorbar(mBetaSquared(i).Year, mBetaSquared(i).mBetaSquared,mBetaSquared(i).Tot,...
    'Marker',mBetaSquared(i).Marker,'MarkerSize',12,'LineWidth',1,'Color','Black',...
    'MarkerEdgeColor',mBetaSquared(i).Color,'MarkerFaceColor',mBetaSquared(i).Color);
    hz(i).MarkerSize = MarkerSize; hz(i).CapSize = 0; %z.FontSize= 20;
    hz(i).LineStyle= 'none';hz(i).LineWidth= 2;legend hide
    xlim([2004 2021])
hold on
end
xticks([2005 2011 mBetaSquared(12).Year-0.5]);
yticks([-3 -2 -1 0 1 2 ]);
ylim([-3.2 2])
set(gca,'XAxisLocation', 'top')
grid off
ylabel(sprintf('{\\it m}_\\nu^2 (eV^2)'));%c^4 
l2=legend([ hz(10) hz(11)  hz(12) pdg2],...
    [mBetaSquared(10).Experiment ' '],...
    [mBetaSquared(11).Experiment ' '],...
    [mBetaSquared(12).Experiment ' '],...
    [mBetaSquaredPDG.Label ' '],...
    'FontSize',LocalFontSize-7,'Location','SouthEast');
l2.NumColumns=2;
legend boxoff;
ylim([-4.5 2.8])
PRLFormat
set(gca,'FontSize',LocalFontSize);
set(gca,'XMinorTick','off');
set(get(gca,'XLabel'),'FontSize',LocalFontSize+2);
set(get(gca,'YLabel'),'FontSize',LocalFontSize+2);
set(gca,'TickDir','out');


if strcmp(SavePlot,'ON')
export_fig(gcf,savefile,'-r300');
export_fig(gcf,strrep(savefile,'.png','.pdf'));
end