function [pHandle,legStr] = KSN1ContourPlot(mnu4Sq,sin2T4,chi2,chi2_ref,CL,varargin)
p = inputParser;
p.addParameter('PlotSplines','ON',@(x) ismember(x,{'ON','OFF'}));
p.addParameter('HoldOn','OFF',@(x) ismember(x,{'ON','OFF'}));
p.addParameter('Color','',@(x) ischar(x) || iscell(x));
p.parse(varargin{:});

PlotSplines = p.Results.PlotSplines;
HoldOn      = p.Results.HoldOn;
Color       = p.Results.Color;

if strcmp(HoldOn,'OFF')
    GetFigure;
end
mnu4Sq_contour = cell(numel(CL),1);
sin2T4_contour = cell(numel(CL),1);
legStr         = cell(numel(CL),1);
for i=1:numel(CL)
    PlotArg = {'-','LineWidth',2};
    if ischar(Color)
        myColor = Color;
    elseif iscell(Color)
        myColor = Color{i};
    end
    if ~isempty(Color)
        PlotArg = {PlotArg{:},'Color',rgb(myColor)};
    end
    [mnu4Sq_contour{i}, sin2T4_contour{i}] = ...
        KSN1Grid2Contour(mnu4Sq,sin2T4,chi2,chi2_ref,CL(i));
    if strcmp(PlotSplines,'OFF')
        pHandle  = plot(sin2T4_contour{i},mnu4Sq_contour{i},...
            PlotArg{:},'MarkerSize',20);
    else
        y = linspace(min(mnu4Sq_contour{i}),max(mnu4Sq_contour{i}),1e3);
        x = interp1(mnu4Sq_contour{i},sin2T4_contour{i},y,'spline');
        pHandle = plot(x,y, PlotArg{:});
    end
    hold on;
    legStr{i} = sprintf('%.0f%% C.L.',CL(i));
end

set(gca,'YScale','log');
set(gca,'XScale','log');

PrettyFigureFormat;
xlabel('|U_{e4}|^2');
ylabel(sprintf('{\\itm}_4 (eV^2)'));

end