% plot script for KSN1 sterile contour plot
% Lisa, April 2020
%% settings
Nathan = 'OFF'; % compare with nathan
CL = 95;
SavePlot = 'ON';
PlotSplines = 'OFF'; % plot smooth spline instead of points
range = 65;%
nGridSteps = 50;
chi2Str = 'chi2CMShape';
DataType = {'Twin','Real'};
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
    'RecomputeFlag','OFF'};


[mnu4Sq_T,sin2T4_T,chi2_T,chi2_ref_T,savefile_T] = KSN1GridSearch(...
    'DataType','Twin',GridArg{:});
[mnu4Sq_D,sin2T4_D,chi2_D,chi2_ref_D,savefile_D] = KSN1GridSearch(...
    'DataType','Real',GridArg{:});
%% plot
[pTwin,legStrTwin] = KSN1ContourPlot('mnu4Sq',mnu4Sq_T,'sin2T4',sin2T4_T,...
    'chi2',chi2_T,'chi2_ref',chi2_ref_T,'CL',CL,...
    'HoldOn','OFF','PlotSplines','ON','Color','DodgerBlue');
hold on;
[pData,legStrData] = KSN1ContourPlot('mnu4Sq',mnu4Sq_D,'sin2T4',sin2T4_D,...
    'chi2',chi2_D,'chi2_ref',chi2_ref_D,'CL',CL,...
      'HoldOn','ON','PlotSplines','ON','Color','Orange');


%% plot contour
if range==40
    xlim([1e-02 0.4])
elseif range>=90
    xlim([1e-03 0.4])
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
plotname = strrep(strrep(strrep(savefile_D,'results','plots'),'.mat','.png'),...
    'Real','RealTwin');

%% nathan contour
if strcmp(Nathan,'ON')
    if range==95
        nrange = 90;
    else
        nrange = range;
    end
    
    nathanfileD = sprintf('%sNathan/NathanContour_%.0feV_%s_%s.mat',savedir,nrange,'Real',chi2Str);
    nathanfileT = sprintf('%sNathan/NathanContour_%.0feV_%s_%s.mat',savedir,nrange,'Twin',chi2Str);
   

    try
        fD = importdata(nathanfileD);
        pN_Data = plot(fD.sith4_X,fD.m4_Y,'-.','LineWidth',2,'Color',rgb('Tomato'));
        
        fT = importdata(nathanfileT);
        pN_Twin = plot(fT.sith4_X,fT.m4_Y,'-.','LineWidth',2,'Color',rgb('Navy'));
      
        plotname = strrep(plotname,'.png','_CompareN.png');
    catch
        fprintf('nathan file %s not available \n',nathanfile)
    end
    
end
%% legend and title
if strcmp(Nathan,'ON')
leg = legend([pN_Data,pData,pN_Twin,pTwin],...
       'Data - Nathan','Data - Lisa',...
       'Twin - Nathan','Twin - Lisa',...
         'EdgeColor',rgb('Silver'),'Location','southwest');
else
    legend([pData,pTwin],...
        sprintf('Data'),...
        sprintf('Twin'),...
        'EdgeColor',rgb('Silver'),'Location','southwest');
end

title(sprintf('%s , %.0f eV at %.0f%% C.L.',chi2Label,range,CL),'FontWeight','normal');

%% save plot
if strcmp(SavePlot,'ON')
    print(plotname,'-dpng','-r450');
    fprintf('save plot to %s \n',plotname);
end
