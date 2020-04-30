% plot script for KSN1 sterile contour plot
% Lisa, April 2020
%% settings
Nathan = 'ON'; % compare with nathan
CL = [90];
SavePlot = 'OFF';
PlotSplines = 'OFF'; % plot smooth spline instead of points
range = 90;%
nGridSteps =50;%50;
chi2Str = 'chi2CMShape';
DataType = 'Real';
freePar = 'E0 Bkg Norm';
RunList = 'KNM1';
SmartGrid = 'OFF';

%% load grid (or calculate if doesn't exist)
[mnu4Sq,sin2T4,chi2,chi2_ref,savefile] = KSN1GridSearch('range',range,...
    'nGridSteps',nGridSteps,...
    'chi2',chi2Str,...
    'DataType',DataType,...
    'freePar',freePar,...
    'RunList',RunList,...
    'SmartGrid',SmartGrid,...
    'RecomputeFlag','OFF');

%% plot contour
[pHandle,legStr] = KSN1ContourPlot(mnu4Sq,sin2T4,chi2,chi2_ref,CL);
if range==40
    xlim([1e-02 0.5])
elseif range>=90
    xlim([1e-03 0.5])
end
ylim([1 (range+5)^2])

if strcmp(chi2Str,'chi2Stat')
    chi2Label = 'stat. only';
else
    chi2Label = 'stat. and syst.';
end

%% plot label
savedir = [getenv('SamakPath'),'ksn1ana/LisaSterile/results/'];
plotdir = strrep(savedir,'results','plots');
MakeDir(plotdir);
plotname = strrep(strrep(savefile,'results','plots'),'.mat','.png');

%% nathan contour
if strcmp(Nathan,'ON')
    if range==95
        nrange = 90;
    else
        nrange = range;
    end
    nathanfile = sprintf('%sNathan/NathanContour_%.0feV_%s_%s.mat',savedir,nrange,DataType,chi2Str);
    try
        d = importdata(nathanfile);
        hold on;
        pN = plot(d.sith4_X,d.m4_Y,'-.','LineWidth',2);
        legStr = {'Lisa','Nathan'};
        plotname = strrep(plotname,'.png','_CompareN.png');
    catch
        fprintf('nathan file %s not available \n',nathanfile)
    end
end

%% legend and title 
 leg = legend(legStr,'EdgeColor',rgb('Silver'),'Location','southwest');
if numel(CL)==1
    title(sprintf('%s , %s , %.0f eV at %.0f%% C.L.',DataType,chi2Label,range,CL),'FontWeight','normal');
else
    title(sprintf('%s , %s , %.0f eV',DataType,chi2Label,range),'FontWeight','normal');
end
%% save plot
if strcmp(SavePlot,'ON')
    print(plotname,'-dpng','-r450');
end
fprintf('save plot to %s \n',plotname);