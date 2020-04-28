function plotHandle = ScatterHist2(x,y,varargin)
% 2-D scatter plot with histograms
p = inputParser;
p.addParameter('xName','',@(x)ischar(x) || isempty(x));
p.addParameter('yName','',@(x)ischar(x) || isempty(x));
p.addParameter('RefLine','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('LegFlag','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('SaveAs','',@(x)ischar(x) || isempty(x));
p.parse(varargin{:});
xName   = p.Results.xName;
yName   = p.Results.yName;
RefLine = p.Results.RefLine;
LegFlag = p.Results.LegFlag;
SaveAs  = p.Results.SaveAs;

% plot
plotHandle = figure('Units','normalized','Position',[0.1,0.1,0.8,0.8]);
scatterhist(x,y,'Direction','out',...);%,'HistogramDisplayStyle','bar',...
    'Location','Northeast','Color',rgb('DodgerBlue'),'Kernel','off');
xlabel(xName);
ylabel(yName);
PrettyFigureFormat('FontSize',26);

if strcmp(RefLine,'ON')
    hold on
    ylimits = ylim;
    xlimits = xlim;
    prefx = plot(mean(x).*ones(10,1),linspace(ylimits(1),ylimits(2),10),'k:','LineWidth',2);
    prefy = plot(linspace(xlimits(1),xlimits(2),10),mean(y).*ones(10,1),'k:','LineWidth',2); 
end

if strcmp(LegFlag,'ON')
    Correlation = corrcoef(x,y);
    leg = legend(sprintf('\\rho = %.3g',Correlation(1,2)));
    leg.EdgeColor = rgb('Silver');
    leg.Location = 'northwest';
end
if ~isempty(SaveAs)
    print(gcf,SaveAs,'-dpng','-r500');
end
end