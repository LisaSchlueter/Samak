
%% Load Data
[mBetaSquared , mBetaSquaredPDG]  = mBetaSquaredExperiments();

%% Plot
LocalFontSize = 18;
MarkerSize = 8;

savefile=sprintf('plots/m2resultsVsyear2019.png');
fig1 = figure('Units','normalized','Position',[0.1,0.1,0.4,1.2]);
s1 = subplot(3,1,1);

[pdg1 pdg2]= boundedline([mBetaSquared(1).Year-1 mBetaSquared(end).Year+1],[-0.6 -0.6],1.9,'alpha','cmap',rgb('CadetBlue'));
hold on
for i=1:12
    hz(i)=errorbar(mBetaSquared(i).Year, mBetaSquared(i).mBetaSquared,mBetaSquared(i).Tot,...
    'Marker',mBetaSquared(i).Marker,'MarkerSize',MarkerSize,'LineWidth',1,'Color','Black',...
    'MarkerEdgeColor',mBetaSquared(i).Color,'MarkerFaceColor',mBetaSquared(i).Color);
    hz(i).CapSize = 0; hz(i).LineStyle= 'none';hz(i).LineWidth= 3;
end


ylabel(sprintf('{\\itm}_\\nu^2 (eV^2)'));%c^4 grid on

PRLFormat
set(gca,'FontSize',LocalFontSize);
set(get(gca,'XLabel'),'FontSize',LocalFontSize+4);
set(get(gca,'YLabel'),'FontSize',LocalFontSize+4);
set(gca,'XMinorTick','off');
xticklabels('')
box on;
ylim([-250 110])
xlim([mBetaSquared(1).Year-1 mBetaSquared(12).Year+1]);
l=legend([hz(10) hz(11)  hz(12) pdg2],...
    [mBetaSquared(10).Experiment ' '],...
    [mBetaSquared(11).Experiment ' '],...
    [mBetaSquared(12).Experiment ' '],...
    [mBetaSquaredPDG.Label ' '],...
    'FontSize',LocalFontSize-7,'Location','northoutside');
l.NumColumns=2;
legend boxoff;

% Zoom
hold on

z = axes(gcf,'position',[0.44 0.73 .45 .1],'FontSize',LocalFontSize-4,'XAxisLocation', 'top');
box on 
[pdg1 pdg2]=boundedline([mBetaSquared(1).Year-1 mBetaSquared(end).Year+1],[-0.6 -0.6],1.9,'alpha','cmap',rgb('CadetBlue'));
hold on;
for i=10:13
    hz(i)=errorbar(mBetaSquared(i).Year, mBetaSquared(i).mBetaSquared,mBetaSquared(i).Tot,...
    'Marker',mBetaSquared(i).Marker,'MarkerSize',MarkerSize,'LineWidth',1,'Color','Black',...
    'MarkerEdgeColor',mBetaSquared(i).Color,'MarkerFaceColor',mBetaSquared(i).Color);
    hz(i).MarkerSize = 16; hz(i).CapSize = 0; %z.FontSize= 20;
    hz(i).LineStyle= 'none';hz(i).LineWidth= 5;legend hide
    xlim([2004 2021])
hold on
end
xticks([2005 2011 mBetaSquared(12).Year-0.5]);
yticks([-3 -2 -1 0 1 2 ]);
ylim([-3.2 2])
set(gca,'XAxisLocation', 'top')
grid off
ylabel(sprintf('{\\it m}_\\nu^2 (eV^2)'));%c^4 

ylim([-4.5 2.8])
PRLFormat
set(gca,'FontSize',LocalFontSize);
set(gca,'XMinorTick','off');
set(get(gca,'XLabel'),'FontSize',LocalFontSize+2);
set(get(gca,'YLabel'),'FontSize',LocalFontSize+2);
% % export_fig(gcf,savefile,'-r300');
% % export_fig(gcf,strrep(savefile,'.png','.pdf'));

%% subplot 2 
% Load Data
[mBetaSquared , mBetaSquaredPDG]  = mBetaSquaredExperiments();

subplot(3,1,2);

stem([mBetaSquared.Year], [mBetaSquared.Sys],'LineStyle','--','Color','Black','LineWidth',2,'Marker','none');
hold on
for i=1:12
    hz(i)=errorbar(mBetaSquared(i).Year, mBetaSquared(i).Sys,mBetaSquared(i).Sys*0,'Marker',mBetaSquared(i).Marker,'MarkerSize',MarkerSize,'LineWidth',1,'Color','Black',...
     'MarkerEdgeColor',mBetaSquared(i).Color,'MarkerFaceColor',mBetaSquared(i).Color);
     hz(i).CapSize = 0; hz(i).LineStyle= 'none';hz(i).LineWidth= 3;
    hold on
end

xlabel('Year');
ylabel(sprintf('{\\itm}_\\nu^2 systematic uncertainty (eV^2)'));

grid on

PRLFormat
set(gca,'FontSize',LocalFontSize);
set(get(gca,'XLabel'),'FontSize',LocalFontSize+4);
set(get(gca,'YLabel'),'FontSize',LocalFontSize+4)
set(gca,'XMinorTick','off');
set(gca,'YScale','log');
box on;
xlim([mBetaSquared(1).Year-1 mBetaSquared(12).Year+1]);
ylim([1e-1 3.3*max([mBetaSquared.Sys])])


%% subplot 3