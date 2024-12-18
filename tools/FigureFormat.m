%% Function that sets default values for matlab text, figure, axis, legend and colorbar objects. As 
% not every property may be accessed through the `set(groot,...)` interface, some are dynamically
% assigned at the end of this function. Hence be sure to call this function at the end of your
% marco.
% The below listed default Property Values were derived from Property Values listed through 
% `get(groot,'factory')`.
%
% Marc Korzeczek, June 2017
% mailto: marc.korzeczek@kit.edu
%

function []=FigureFormat(varargin)
p=inputParser();
p.addParameter('FontWeight','normal',@(x)ismember(x,{'normal','bold'}));
p.addParameter('FontSize',12,@(x)isfloat(x)&&x>0);
p.addParameter('FigurePosition',[5 5 15 10],@(x)isfloat(x));
p.addParameter('TextInterpreter','latex',@(x)ismember(x,{'tex','latex'}));
p.addParameter('MajorGrid','on',@(x)ismember(x,{'on','off'}));
p.addParameter('MinorGrid','on',@(x)ismember(x,{'on','off'}));
p.parse(varargin{:});

TextInterpreter = p.Results.TextInterpreter;
FontWeight  = p.Results.FontWeight;
FontSize  = p.Results.FontSize;
FigurePosition = p.Results.FigurePosition;
MajorGrid = p.Results.MajorGrid;
MinorGrid = p.Results.MinorGrid;


%% Define default text properties
TextProperties={ ...
    'defaultTextBackgroundColor', 'none',...
    'defaultTextBusyAction', 'queue',...
    'defaultTextButtonDownFcn', '',...
    'defaultTextClipping', 'off',...
    'defaultTextColor', [0 0 0],...
    'defaultTextCreateFcn', '',...
    'defaultTextDeleteFcn', '',...
    'defaultTextEdgeColor', 'none',...
    'defaultTextEditing', 'off',...
    'defaultTextFontAngle', 'normal',...
    'defaultTextFontName', 'Helvetica',...
    'defaultTextFontSize', FontSize,...
    'defaultTextFontSmoothing', 'on',...
    'defaultTextFontUnits', 'points',...
    'defaultTextFontWeight', 'normal',...
    'defaultTextHandleVisibility', 'on',...
    'defaultTextHitTest', 'on',...
    'defaultTextHorizontalAlignment', 'left',...
    'defaultTextInterpreter', TextInterpreter,...
    'defaultTextInterruptible', 'on',...
    'defaultTextLineStyle', '-',...
    'defaultTextLineWidth', 0.5000,...
    'defaultTextMargin', 3,...
    ...%'defaultTextParent', [0�0 GraphicsPlaceholder],...
    'defaultTextPickableParts', 'visible',...
    ...%'defaultTextPosition', [3�1 double],...
    'defaultTextRotation', 0,...
    'defaultTextSelected', 'off',...
    'defaultTextSelectionHighlight', 'on',...
    'defaultTextString', '',...
    'defaultTextTag', '',...
    'defaultTextUnits', 'data',...
    'defaultTextUserData', [],...
    'defaultTextVerticalAlignment', 'middle',...
    'defaultTextVisible', 'on'};


%% Define default figure properties
FigureProperties={ ...
    ...%''defaultFigureAlphamap', [1�64 double],...
    'defaultFigureBusyAction', 'queue',...
    'defaultFigureButtonDownFcn', '',...
    'defaultFigureClipping', 'on',...
    'defaultFigureCloseRequestFcn', 'closereq',...
    'defaultFigureColor', [0.9400 0.9400 0.9400],...
    ...%''defaultFigureColormap', [64�3 double],...
    'defaultFigureCreateFcn', '',...
    'defaultFigureDeleteFcn', '',...
    'defaultFigureDockControls', 'on',...
    'defaultFigureFileName', '',...
    'defaultFigureGraphicsSmoothing', 'on',...
    'defaultFigureHandleVisibility', 'on',...
    'defaultFigureIntegerHandle', 'on',...
    'defaultFigureInterruptible', 'on',...
    'defaultFigureInvertHardcopy', 'on',...
    'defaultFigureKeyPressFcn', '',...
    'defaultFigureKeyReleaseFcn', '',...
    'defaultFigureMenuBar', 'figure',...
    'defaultFigureName', '',...
    'defaultFigureNextPlot', 'add',...
    'defaultFigureNumberTitle', 'on',...
    'defaultFigurePaperOrientation', 'portrait',...
    'defaultFigurePaperPosition', FigurePosition,...
    'defaultFigurePaperPositionMode', 'auto',...
    'defaultFigurePaperSize', [FigurePosition(3) FigurePosition(4)]*1.05,...
    'defaultFigurePaperType', 'usletter',...
    'defaultFigurePaperUnits', 'centimeters',...
    'defaultFigurePointer', 'arrow',...
    ...%'defaultFigurePointerShapeCData', [16�16 double],...
    'defaultFigurePointerShapeHotSpot', [1 1],...
    'defaultFigurePosition', FigurePosition*1.05,...
    'defaultFigureRenderer', 'opengl',...
    'defaultFigureRendererMode', 'auto',...
    'defaultFigureResize', 'on',...
    'defaultFigureSizeChangedFcn', '',...
    'defaultFigureTag', '',...
    'defaultFigureToolBar', 'auto',...
    'defaultFigureUnits', 'centimeters',...
    'defaultFigureUserData', [],...
    'defaultFigureVisible', 'on',...
    'defaultFigureWindowButtonDownFcn', '',...
    'defaultFigureWindowButtonMotionFcn', '',...
    'defaultFigureWindowButtonUpFcn', '',...
    'defaultFigureWindowKeyPressFcn', '',...
    'defaultFigureWindowKeyReleaseFcn', '',...
    'defaultFigureWindowScrollWheelFcn', '',...
    'defaultFigureWindowStyle', 'normal'};


%% Define default axes properties
AxesProperties={ ...
    'defaultAxesALim', [0 1],...
    'defaultAxesALimMode', 'auto',...
    'defaultAxesActivePositionProperty', 'outerposition',...
    'defaultAxesAmbientLightColor', [1 1 1],...
    'defaultAxesBox', 'off',...
    'defaultAxesBoxStyle', 'back',...
    'defaultAxesBusyAction', 'queue',...
    'defaultAxesButtonDownFcn', '',...
    'defaultAxesCLim', [0 1],...
    'defaultAxesCLimMode', 'auto',...
    'defaultAxesCameraPosition', [0 0 1],...
    'defaultAxesCameraPositionMode', 'auto',...
    'defaultAxesCameraTarget', [0 0 0],...
    'defaultAxesCameraTargetMode', 'auto',...
    'defaultAxesCameraUpVector', [0 1 0],...
    'defaultAxesCameraUpVectorMode', 'auto',...
    'defaultAxesCameraViewAngle', 6.6086,...
    'defaultAxesCameraViewAngleMode', 'auto',...
    'defaultAxesClipping', 'on',...
    'defaultAxesClippingStyle', '3dbox',...
    'defaultAxesColor', [1 1 1],...
    ...%'defaultAxesColorOrder', [7�3 double],...
    'defaultAxesColorOrderIndex', 1,...
    'defaultAxesCreateFcn', '',...
    'defaultAxesDataAspectRatio', [1 1 1],...
    'defaultAxesDataAspectRatioMode', 'auto',...
    'defaultAxesDeleteFcn', '',...
    'defaultAxesFontAngle', 'normal',...
    'defaultAxesFontName', 'Helvetica',...
    'defaultAxesFontSize', FontSize,...
    'defaultAxesFontSmoothing', 'on',...
    'defaultAxesFontUnits', 'points',...
    'defaultAxesFontWeight', 'normal',...
    'defaultAxesGridAlpha', 0.2000,...
    'defaultAxesGridAlphaMode', 'auto',...
    'defaultAxesGridColor', [0.1500 0.1500 0.1500],...
    'defaultAxesGridColorMode', 'auto',...
    'defaultAxesGridLineStyle', '-',...
    'defaultAxesHandleVisibility', 'on',...
    'defaultAxesHitTest', 'on',...
    'defaultAxesInterruptible', 'on',...
    'defaultAxesLabelFontSizeMultiplier', 1.1000,...
    'defaultAxesLayer', 'bottom',...
    ...%'defaultAxesLegend', [0�0 GraphicsPlaceholder],...
    'defaultAxesLineStyleOrder', '-',...
    'defaultAxesLineStyleOrderIndex', 1,...
    'defaultAxesLineWidth', 0.5000,...
    'defaultAxesMinorGridAlpha', 0.2500,...
    'defaultAxesMinorGridAlphaMode', 'auto',...
    'defaultAxesMinorGridColor', [0.1000 0.1000 0.1000],...
    'defaultAxesMinorGridColorMode', 'auto',...
    'defaultAxesMinorGridLineStyle', ':',...
    'defaultAxesNextPlot', 'replace',...
    'defaultAxesOuterPosition', [0 0 1 1],...
    ...%'defaultAxesParent', [0�0 GraphicsPlaceholder],...
    'defaultAxesPickableParts', 'visible',...
    'defaultAxesPlotBoxAspectRatio', [1 1 1],...
    'defaultAxesPlotBoxAspectRatioMode', 'auto',...
    'defaultAxesPosition', [0.1300 0.1100 0.7750 0.8150],...
    'defaultAxesProjection', 'orthographic',...
    'defaultAxesSelected', 'off',...
    'defaultAxesSelectionHighlight', 'on',...
    'defaultAxesSortMethod', 'depth',...
    'defaultAxesTag', '',...
    'defaultAxesTickDir', 'in',...
    'defaultAxesTickDirMode', 'auto',...
    'defaultAxesTickLabelInterpreter', TextInterpreter,...
    'defaultAxesTickLength', [0.0100 0.0250],...
    ...%'defaultAxesTitle', [0�0 Text],...
    'defaultAxesTitleFontSizeMultiplier', 1.1000,...
    'defaultAxesTitleFontWeight', 'bold',...
    'defaultAxesUnits', 'normalized',...
    'defaultAxesUserData', [],...
    'defaultAxesView', [0 90],...
    'defaultAxesVisible', 'on',...
    ...%'defaultAxesXAxis', [0�0 GraphicsPlaceholder],...
    'defaultAxesXAxisLocation', 'bottom',...
    'defaultAxesXColor', [0.1500 0.1500 0.1500],...
    'defaultAxesXColorMode', 'auto',...
    'defaultAxesXDir', 'normal',...
    'defaultAxesXGrid', MajorGrid,...
    ...%'defaultAxesXLabel', [0�0 GraphicsPlaceholder],...
    'defaultAxesXLim', [0 1],...
    'defaultAxesXLimMode', 'auto',...
    'defaultAxesXMinorGrid', MinorGrid,...
    'defaultAxesXMinorTick', 'off',...
    'defaultAxesXScale', 'linear',...
    'defaultAxesXTick', [],...
    'defaultAxesXTickLabel', '',...
    'defaultAxesXTickLabelMode', 'auto',...
    'defaultAxesXTickLabelRotation', 0,...
    'defaultAxesXTickMode', 'auto',...
    'defaultAxesYAxisLocation', 'left',...
    'defaultAxesYColor', [0.1500 0.1500 0.1500],...
    'defaultAxesYColorMode', 'auto',...
    'defaultAxesYDir', 'normal',...
    'defaultAxesYGrid', MajorGrid,...
    ...%'defaultAxesYLabel', [0�0 GraphicsPlaceholder],...
    'defaultAxesYLim', [0 1],...
    'defaultAxesYLimMode', 'auto',...
    'defaultAxesYMinorGrid', MinorGrid,...
    'defaultAxesYMinorTick', 'off',...
    'defaultAxesYScale', 'linear',...
    'defaultAxesYTick', [],...
    'defaultAxesYTickLabel', '',...
    'defaultAxesYTickLabelMode', 'auto',...
    'defaultAxesYTickLabelRotation', 0,...
    'defaultAxesYTickMode', 'auto',...
    ...%'defaultAxesZAxis', [0�0 GraphicsPlaceholder],...
    'defaultAxesZColor', [0.1500 0.1500 0.1500],...
    'defaultAxesZColorMode', 'auto',...
    'defaultAxesZDir', 'normal',...
    'defaultAxesZGrid', MajorGrid,...
    ...%'defaultAxesZLabel', [0�0 GraphicsPlaceholder],...
    'defaultAxesZLim', [0 1],...
    'defaultAxesZLimMode', 'auto',...
    'defaultAxesZMinorGrid', MinorGrid,...
    'defaultAxesZMinorTick', 'off',...
    'defaultAxesZScale', 'linear',...
    'defaultAxesZTick', [],...
    'defaultAxesZTickLabel', '',...
    'defaultAxesZTickLabelMode', 'auto',...
    'defaultAxesZTickLabelRotation', 0,...
    'defaultAxesZTickMode', 'auto'};


%% Define default legend properties
LegendProperties={ ...
    'defaultLegendAutoUpdate', 'on',...
    'defaultLegendBox', 'on',...
    'defaultLegendBusyAction', 'queue',...
    'defaultLegendButtonDownFcn', '',...
    'defaultLegendColor', [1 1 1],...
    'defaultLegendCreateFcn', '',...
    'defaultLegendDeleteFcn', '',...
    'defaultLegendEdgeColor', [0.2000 0.2000 0.2000],...
    'defaultLegendFontAngle', 'normal',...
    'defaultLegendFontName', 'Helvetica',...
    'defaultLegendFontSize', FontSize-1,...
    'defaultLegendFontWeight', 'normal',...
    'defaultLegendHandleVisibility', 'on',...
    'defaultLegendHitTest', 'on',...
    'defaultLegendInterpreter', TextInterpreter,...
    'defaultLegendInterruptible', 'on',...
    'defaultLegendItemHitFcn', @defaultItemHitCallback,...
    'defaultLegendLineWidth', 0.5000,...
    'defaultLegendLocation', 'northeast',...
    'defaultLegendOrientation', 'vertical',...
    ...%'defaultLegendParent', [0�0 GraphicsPlaceholder],...
    'defaultLegendPickableParts', 'visible',...
    'defaultLegendPosition', [0 0 1 1],...
    'defaultLegendSelected', 'off',...
    'defaultLegendSelectionHighlight', 'on',...
    'defaultLegendString', [],...
    'defaultLegendTag', '',...
    'defaultLegendTextColor', [0 0 0],...
    ...%'defaultLegendTitle', [0�0 Text],...
    'defaultLegendUnits', 'normalized',...
    'defaultLegendUserData', [],...
    'defaultLegendVisible', 'on'};

%% Define default colorbar properties
ColorbarProperties={ ...
    'defaultColorbarAxisLocation', 'out' ,...
    'defaultColorbarAxisLocationMode', 'auto',...
    'defaultColorbarBox', 'on',...
    'defaultColorbarBusyAction', 'queue',...
    'defaultColorbarButtonDownFcn', '',...
    'defaultColorbarColor', [0 0 0],...
    'defaultColorbarCreateFcn', '',...
    'defaultColorbarDeleteFcn', '',...
    'defaultColorbarDirection', 'normal',...
    'defaultColorbarFontAngle', 'normal',...
    'defaultColorbarFontName', 'Helvetica',...
    'defaultColorbarFontSize', FontSize-1,...
    'defaultColorbarFontWeight', 'normal',...
    'defaultColorbarHandleVisibility', 'on',...
    'defaultColorbarHitTest', 'on',...
    'defaultColorbarInterruptible', 'on',...
    ...%'defaultColorbarLabel', [0�0 matlab.mixin.Heterogeneous],...
    'defaultColorbarLimits', [0 1],...
    'defaultColorbarLimitsMode', 'auto',...
    'defaultColorbarLineWidth', 0.5000,...
    'defaultColorbarLocation', 'eastoutside',...
    ...%'defaultColorbarParent', [0�0 GraphicsPlaceholder],...
    'defaultColorbarPickableParts', 'visible',...
    'defaultColorbarPosition', [0 0 1 1],...
    'defaultColorbarSelected', 'off',...
    'defaultColorbarSelectionHighlight', 'on',...
    'defaultColorbarTag', '',...
    'defaultColorbarTickDirection', 'in',...
    'defaultColorbarTickLabelInterpreter', TextInterpreter,...
    'defaultColorbarTickLabels', '',...
    'defaultColorbarTickLabelsMode', 'auto',...
    'defaultColorbarTickLength', 0.0100,...
    'defaultColorbarTicks', [],...
    'defaultColorbarTicksMode', 'auto',...
    'defaultColorbarUnits', 'normalized',...
    'defaultColorbarUserData', [],...
    'defaultColorbarVisible', 'on'};



%% Statically set default properties
set(groot,TextProperties{:}, FigureProperties{:}, AxesProperties{:}, ...
    LegendProperties{:}, ColorbarProperties{:});



%% Dynamically set properties
figHandles = findall(0,'Type','figure');
for i=1:numel(figHandles)
    % modify figure properties
    fig=figHandles(i);
    % modify axis properties
    allaxes = findall(fig, 'type', 'axes');
    for j=1:numel(allaxes)
        ax=allaxes(j);
        set(ax,'LooseInset',get(ax,'TightInset'));
    end
    
    % modify legend properties
    %alllegend = findall(fig, 'type', 'legend');
    %for j=1:numel(alllegend)
        %leg=alllegend(j);
        %set(leg,'Property',Value);   
    %end
    
    % modify colorbar properties
    allcolorbar = findall(fig, 'type', 'colorbar');
    for j=1:numel(allcolorbar)
        col=allcolorbar(j);
        set(get(col,'Label'),'FontSize',FontSize,'FontWeight',FontWeight,...
            'Interpreter',TextInterpreter);
    end
    
    refresh(fig);
end


end