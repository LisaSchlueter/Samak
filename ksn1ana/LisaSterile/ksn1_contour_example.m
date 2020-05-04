% example script 

%% settings
CL = 0.95;
range = 95;%
nGridSteps = 50;
chi2Str = 'chi2CMShape';
DataType = 'Real';
freePar = 'mNu E0 Bkg Norm';
RunList = 'KNM1';
SmartGrid = 'OFF';
pullFlag = 12;
%% plot options
PlotContour = 'OFF';
PlotGrid    = 'ON';
if strcmp(chi2Str,'chi2Stat')
    chi2Label = 'stat. only';
else
    chi2Label = 'stat. and syst.';
end
if strcmp(DataType,'Real')
    DataLabel = 'Data';
else
    DataLabel = 'Twins';
end
titleStr = sprintf('%s (%s) %.0f eV range',DataLabel,chi2Label,range);
%% load grid (or calculate if doesn't exist)
[mnu4Sq,sin2T4,chi2,chi2_ref,savefile] = KSN1GridSearch('range',range,...
    'nGridSteps',nGridSteps,...
    'chi2',chi2Str,...
    'DataType',DataType,...
    'freePar',freePar,...
    'RunList',RunList,...
    'SmartGrid',SmartGrid,...
    'RecomputeFlag','OFF',...
    'pullFlag',pullFlag);

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
plotnameGrid = strrep(strrep(extractAfter(savefile,'results/'),'.mat','.png'),'KSN1','Chi2Map_KSN1');
if strcmp(PlotGrid,'ON')
     KSN1GridPlot(PlotArg{:},'nInter',1e3,...
    'ContourPlot','Fitrium','SaveAs',plotnameGrid);
end
