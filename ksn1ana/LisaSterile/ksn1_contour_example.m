% example script 

%% settings
range = 95;%
nGridSteps = 50;
chi2Str = 'chi2CMShape';
DataType = 'Real';
freePar = 'E0 Bkg Norm';
RunList = 'KNM1';
SmartGrid = 'OFF';

%% plot options
PlotContour = 'OFF';
PlotGrid    = 'ON';
CL = 0.82;
titleStr = sprintf('%s (%s) %.0f eV range',DataType,chi2Str,range);
%% load grid (or calculate if doesn't exist)
[mnu4Sq,sin2T4,chi2,chi2_ref,savefile] = KSN1GridSearch('range',range,...
    'nGridSteps',nGridSteps,...
    'chi2',chi2Str,...
    'DataType',DataType,...
    'freePar','E0 Bkg Norm',...
    'RunList',RunList,...
    'SmartGrid',SmartGrid,...
    'RecomputeFlag','OFF');

%% plot contour
 PlotArg ={'mnu4Sq',mnu4Sq,...
        'sin2T4',sin2T4,...
        'chi2',chi2,'chi2_ref',chi2_ref,...
        'CL',CL,...
        'titleStr',titleStr};
if strcmp(PlotContour,'ON')
    KSN1ContourPlot(PlotArg{:},'LineStyle','-','PlotSplines','ON');
end

%% plot grid
if strcmp(PlotGrid,'ON')
     KSN1GridPlot(PlotArg{:});
end
