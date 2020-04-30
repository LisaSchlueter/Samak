function [pHandle,legStr] = KSN1ContourPlot(mnu4Sq,sin2T4,chi2,chi2_ref,CL)
PlotSplines = 'OFF';
GetFigure;
mnu4Sq_contour = cell(numel(CL),1);
sin2T4_contour = cell(numel(CL),1);
legStr         = cell(numel(CL),1);
for i=1:numel(CL)
    [mnu4Sq_contour{i}, sin2T4_contour{i}] = ...
        KSN1Grid2Contour(mnu4Sq,sin2T4,chi2,chi2_ref,CL(i));
    if strcmp(PlotSplines,'OFF')
        pHandle  = plot(sin2T4_contour{i},mnu4Sq_contour{i},'-','LineWidth',2,'MarkerSize',20);
    else
        y = linspace(min(mnu4Sq_contour{i}),max(mnu4Sq_contour{i}),1e3);
        x = interp1(mnu4Sq_contour{i},sin2T4_contour{i},y,'spline');
        pHandle = plot(x,y,'-','LineWidth',2);
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