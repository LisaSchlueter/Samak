function [pHandle,legStr] = KSN1ContourPlot(varargin)
p = inputParser;
p.addParameter('mnu4Sq','',@(x)isfloat(x));
p.addParameter('sin2T4','',@(x)isfloat(x));
p.addParameter('chi2','',@(x)isfloat(x));
p.addParameter('chi2_ref','',@(x)isfloat(x));
p.addParameter('CL',0.9,@(x)isfloat(x));
p.addParameter('PlotSplines','ON',@(x) ismember(x,{'ON','OFF'}));
p.addParameter('HoldOn','OFF',@(x) ismember(x,{'ON','OFF'}));
p.addParameter('Color','',@(x) ischar(x) || iscell(x));
p.addParameter('LineStyle','-',@(x) ischar(x));
p.addParameter('titleStr','',@(x)ischar(x));
p.addParameter('nInter',1e3,@(x)isfloat(x));
p.addParameter('BestFit','OFF',@(x) ismember(x,{'ON','OFF'}));
p.addParameter('Method','New',@(x) ismember(x,{'New','Old'}));

p.parse(varargin{:});

mnu4Sq      = p.Results.mnu4Sq;
sin2T4      = p.Results.sin2T4;
chi2        = p.Results.chi2;
chi2_ref    = p.Results.chi2_ref;
CL          = p.Results.CL;
PlotSplines = p.Results.PlotSplines;
HoldOn      = p.Results.HoldOn;
Color       = p.Results.Color;
LineStyle   = p.Results.LineStyle;
titleStr    = p.Results.titleStr;
nInter      = p.Results.nInter;
BestFit     = p.Results.BestFit;
Method      = p.Results.Method;
if isempty(mnu4Sq)
    load('KSN1_SmartGridInit_40eVrange.mat','mnu4Sq','sin2T4','chi2','chi2_ref');
    fprintf('no input grid specified - load default grid \n')
    titleStr = 'Default grid 40 eV range (twins)';
end
%% plot
if strcmp(HoldOn,'OFF')
    GetFigure;
end
mnu4Sq_contour = cell(numel(CL),1);
sin2T4_contour = cell(numel(CL),1);
legStr         = cell(numel(CL),1);
for i=1:numel(CL)
    PlotArg = {'LineWidth',2,'MarkerSize',4};
    if ischar(Color)
        myColor = Color;
    elseif iscell(Color)
        myColor = Color{i};
    end
    if ~isempty(Color)
        PlotArg = [PlotArg,{'Color',rgb(myColor)}];
    end
 
    [mnu4Sq_contour{i}, sin2T4_contour{i}] = ...
        KSN1Grid2Contour(mnu4Sq,sin2T4,chi2,chi2_ref,CL(i),'Mode',Method);
    
    if strcmp(Method,'New')
        mnu4Sq_tmp = mnu4Sq_contour{i};
        sin2T4_tmp = sin2T4_contour{i};
    else
        mnu4Sq_tmp = mnu4Sq_contour(i);
        sin2T4_tmp = sin2T4_contour(i);
    end
    nContour = size(mnu4Sq_tmp,1);
    if nContour>1 % || numel(mnu4Sq_tmp{i})>1e4
        PlotSplines='OFF';
    end
    
    for c=1:nContour
        if strcmp(PlotSplines,'OFF')
            pHandle  = plot(sin2T4_tmp{c}, mnu4Sq_tmp{c},...
                LineStyle,PlotArg{:},'MarkerSize',10);
        else
            try 
                y = logspace(log10(min(mnu4Sq_tmp{c})),log10(max(mnu4Sq_tmp{c})),1e3);
                x = interp1(mnu4Sq_tmp{c},sin2T4_tmp{c},y,'spline');
                pHandle = plot(x,y,LineStyle,PlotArg{:});
            catch
                fprintf('Interpolation failed/interupted - probably closed contour \n')
                pHandle  = plot(sin2T4_tmp{c}, mnu4Sq_tmp{c},...
                    LineStyle,PlotArg{:},'MarkerSize',10);
            end
        end
        hold on;
    end
    legStr{i} = sprintf('%.0f%% C.L.',CL(i));
    
end

PrettyFigureFormat;

%% plot best fit
if strcmp(BestFit,'ON')
    
    [row, col] = find(chi2 == min(chi2(:)));
    mnu4Sq_bf =  mnu4Sq(col,row);
    sin2T4_bf = sin2T4(col,row);
    
    hold on;
    plot(sin2T4_bf,mnu4Sq_bf,'x','MarkerSize',10);%'Color',rgb(myColor)
end
title(titleStr,'FontWeight','normal','FontSize',get(gca,'FontSize'))

%% axis
set(gca,'YScale','log');
set(gca,'XScale','log');

xlabel('|U_{e4}|^2');
ylabel(sprintf('{\\itm}_4^2 (eV^2)'));

rangeApprox = sqrt(max(max(mnu4Sq)));
if rangeApprox==40
    xlim([1e-02 0.5])
else
    xlim([1e-03 0.5])
end
ylim([1 (rangeApprox)^2]);

end