% plot script for KSN1 sterile contour plot
% Lisa, April 2020
%% settings
Nathan = 'OFF'; % compare with nathan
Fitrium = 'OFF';
Giunti = 'ON';
CL = [95];%, 95 99];
SavePlot = 'ON';
PlotSplines = 'OFF'; % plot smooth spline instead of points
range = 40;%
nGridSteps = 50;
chi2Str = 'chi2Stat';%CMShape';
DataType = 'Twin';
freePar = 'E0 Bkg Norm';
RunList = 'KNM1';
SmartGrid = 'OFF';

%% load grid (or calculate if doesn't exist)
GridArg = {'range',range,...
    'nGridSteps',nGridSteps,...
    'chi2',chi2Str,...
    'freePar',freePar,...
    'RunList',RunList,...
    'SmartGrid',SmartGrid,...
    'RecomputeFlag','OFF',...
     };

[mnu4Sq,sin2T4,chi2,chi2_ref,savefile] = KSN1GridSearch(...
    'DataType',DataType,GridArg{:});
%% collect plot argument

if strcmp(chi2Str,'chi2Stat')
    chi2Label = 'stat. only';
else
    chi2Label = 'stat. and syst.';
end
if numel(CL)==1
    titleStr = sprintf('%s , %s , %.0f eV at %.0f%% C.L.',DataType,chi2Label,range,CL);
else
    titleStr = sprintf('%s , %s , %.0f eV',DataType,chi2Label,range);
end

PlotArg ={'mnu4Sq',mnu4Sq,...
    'sin2T4',sin2T4,...
    'chi2',chi2,'chi2_ref',chi2_ref,...
    'CL',CL,...
    'titleStr',titleStr};
%% plot
[pHandle,legStr] = KSN1ContourPlot(PlotArg{:},'PlotSplines','ON');



%% plot label
savedir = [getenv('SamakPath'),'ksn1ana/LisaSterile/results/'];
plotdir = strrep(savedir,'results','plots');
MakeDir(plotdir);
plotname = strrep(strrep(savefile,'results','plots'),'.mat','.png');

%% range
    if range==95
        nrange = 90;
    else
        nrange = range;
    end
%% nathan
    if strcmp(Nathan,'ON')
        nathanfile = sprintf('%sNathan/NathanContour_%.0feV_%s_%s.mat',savedir,nrange,DataType,chi2Str);
        legStr = {'Lisa','Nathan'};
        extraStr = '_N';
        try
        d = importdata(nathanfile);
        hold on;
        pN = plot(d.sith4_X,d.m4_Y,'-.','LineWidth',2);
        plotname = strrep(plotname,'.png',sprintf('_Compare%s.png',extraStr));
    catch
        fprintf('nathan file %s not available \n',nathanfile)
        end
    end
    
    %% fitrum
    if strcmp(Fitrium,'ON')
    else
        fitriumfile = sprintf('%sOthers/coord_fitrium_40_Twin_95_stat.mat',savedir);
        legStr = {'Samak (Lisa)','Fitrium'};
        extraStr = '_F';
    end
    
    try
        d = importdata(fitriumfile);
        hold on;
        pF = plot(d.sith4_X,d.m4_Y,'-.','LineWidth',2);
        plotname = strrep(plotname,'.png',sprintf('_Compare%s.png',extraStr));
    catch
        fprintf('fitrium file %s not available \n',nathanfile)
    end

%% giunti
    if strcmp(Giunti,'ON')
    else
        giuntifile = sprintf('%sOthers/coord_Giunti.mat',savedir);
        legStr = {'Samak (Lisa)','Fitrium'};
        extraStr = '_G';
    end
    
    try
        d = importdata(giuntifile);
        hold on;
        pF = plot(d.sith4_X,d.m4_Y,'-.','LineWidth',2);
        plotname = strrep(plotname,'.png',sprintf('_Compare%s.png',extraStr));
    catch
        fprintf('giunti file %s not available \n',nathanfile)
    end


%% legend
leg = legend(legStr,'EdgeColor',rgb('Silver'),'Location','southwest');
%% save plot
if strcmp(SavePlot,'ON')
    print(plotname,'-dpng','-r450');
end
fprintf('save plot to %s \n',plotname);