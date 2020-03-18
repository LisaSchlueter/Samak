function E0 = knm2FS_GetE0Twins(varargin)
% define runwise endpoint for KNM2 figure skating twins
% March 2020, Lisa
p = inputParser;
p.addParameter('SanityPlot','ON',@(x)ismember(x,{'ON','OFF'}));
p.parse(varargin{:});
SanityPlot  = p.Results.SanityPlot;

RunListDir = [getenv('SamakPath'),'inputs/RunLists/KNM2/'];    
runsRW1 = sort(importdata([RunListDir,'KNM2_RW1.dat']));
runsRW2 = sort(importdata([RunListDir,'KNM2_RW2.dat']));
runsRW3 = sort(importdata([RunListDir,'KNM2_RW3.dat']));

E0 = zeros(numel([runsRW1,runsRW2,runsRW3]),1);

E0(1:numel(runsRW1)) = 18573.62;
E0(numel(runsRW1)+1:numel(runsRW2)+numel(runsRW1)) = 18573.78;
E0(numel(runsRW2)+numel(runsRW1)+1:end) = 18573.71;


if strcmp(SanityPlot,'ON')
    fE0 = figure('Units','normalized','Position',[0.1,0.1,0.65,0.55]);
    s1 = subplot(1,8,1:6);
    plot(1:numel(E0),E0-mean(E0),'s','MarkerSize',3,'MarkerFaceColor',rgb('DodgerBlue'),'Color',rgb('DodgerBlue'));
    PrettyFigureFormat('FontSize',24);
    xlabel('Run number');
    ylabel(sprintf('E_0^{MC} - \\langleE_0^{MC}\\rangle  (eV)'));
    xlim([0 numel(E0)+1]);
    ax = gca;
    ax.Position = [ax.Position(1) ax.Position(2) ax.Position(3)+0.023 ax.Position(4)];
    
    leg = legend(sprintf(' \\langleE_0^{MC}\\rangle = %.2f eV \n\\sigma(E_0^{MC}) = %.0f meV',mean(E0),1e3.*std(E0)),...
        'Location','northwest');
    leg.EdgeColor = rgb('Silver');
    
    s2 = subplot(1,8,7:8);
    BinWidth = 0.01;
    h1 = histogram(E0-mean(E0));
    h1.Orientation='horizontal';
    h1.BinWidth = BinWidth;
    h1.FaceColor = rgb('DodgerBlue');
    h1.FaceAlpha = 1;
    PrettyFigureFormat;
    
    % get rid of box and x axis
    box off
    set(get(gca,'XAxis'),'Visible','off')
    set(get(gca,'YAxis'),'Visible','off')
    
    linkaxes([s1,s2],'y');
    
    savedir  = [getenv('SamakPath'),'knm2ana/knm2_FigureSkating/plots/'];
    savename = sprintf('%sKnm2TwinsE0.pdf',savedir);
    export_fig(fE0,savename);
end
end