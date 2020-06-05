function m2statVsyear2019(varargin)
p=inputParser;
p.addParameter('Combi','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('SavePlot','OFF',@(x)ismember(x,{'ON','OFF'}));

p.parse(varargin{:});

Combi    = p.Results.Combi;
SavePlot = p.Results.SavePlot;
% Plot of Squared neutrino mass stat/syst errors obtained from tritium beta decay
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
maintitle=sprintf('Squared neutrino mass Uncertainties obtained from tritium \\beta -decay in the period 1990-2019');
savefile=sprintf('plots/m2statVsyear2019.png');
if strcmp(Combi,'OFF')
fig1 = figure('Name','Test','NumberTitle','off','rend','painters','pos',[10 10 1400 600]);
end
% a=annotation('textbox', [0 0.9 1 0.1], ...
%     'String', maintitle, ...
%     'EdgeColor', 'none', ...
%     'HorizontalAlignment', 'center');
% a.FontSize=26;a.FontWeight='normal';
stem([mBetaSquared.Year], [mBetaSquared.Stat],'LineStyle',':','Color','Black','LineWidth',2,'Marker','none');
hold on
for i=1:12
    h(i)=errorbar(mBetaSquared(i).Year, mBetaSquared(i).Stat,mBetaSquared(i).Stat*0,'Marker',mBetaSquared(i).Marker,'MarkerSize',MarkerSize,'LineWidth',1,'Color','Black',...
     'MarkerEdgeColor',mBetaSquared(i).Color,'MarkerFaceColor',mBetaSquared(i).Color);
     h(i).CapSize = 0; h(i).LineStyle= 'none';h(i).LineWidth= 3;
    hold on
end
if strcmp(Combi,'OFF')
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
    'FontSize',13.5,'Location','north');

l.NumColumns=4;
set(l.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.8]));
l.EdgeColor = 'none';
end
xlabel('Year');
ylabel(sprintf('{\\itm}_\\nu^2 statistical uncertainty (eV^2)'));

grid on
% %%FontName = 'Arial';
% LocalFontSize = 22;
% set(gca,'FontName',FontName,'FontSize',LocalFontSize);
% set(get(gca,'XLabel'),'FontName',FontName,'FontSize',LocalFontSize);
% set(get(gca,'YLabel'),'FontName',FontName,'FontSize',LocalFontSize);
% set(gca,'TickLength',[.01 .01]);
% set(gca,'XMinorTick','off');
% set(gca,'YMinorTick','on');
% set(gcf,'Color','w');
% set(gca,'YScale','log');

PRLFormat
set(gca,'FontSize',LocalFontSize);
set(get(gca,'XLabel'),'FontSize',LocalFontSize+4);
set(get(gca,'YLabel'),'FontSize',LocalFontSize+4)
set(gca,'XMinorTick','off');
set(gca,'YScale','log');
box on;
xlim([mBetaSquared(1).Year-1 mBetaSquared(12).Year+1]);
ylim([1e-1 3.3*max([mBetaSquared.Stat])])
if strcmp(SavePlot,'ON')
export_fig(gcf,savefile,'-r300');
export_fig(gcf,strrep(savefile,'.png','.pdf'));
end