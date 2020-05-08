% contour with randomized twins
%% settings
CL = 0.95;
range = 95;%
nGridSteps = 25;
chi2Str = 'chi2CMShape';
DataType = 'Twin';
freePar = 'E0 Bkg Norm';
RunList = 'KNM1';
SmartGrid = 'OFF';
RandMC = [1:500];%,500:643];
SysBudget = 24;
%% init
nContours = numel(RandMC);
mnu4Sq   = cell(nContours,1);
sin2T4   = cell(nContours,1);
chi2     = cell(nContours,1);
chi2_ref = cell(nContours,1);
savefile = cell(nContours,1);
FitResults_ref = cell(nContours,1);
DeltaChi2 = zeros(nContours,1);
%% load grid (or calculate if doesn't exist)
mnu4Sq_bf = zeros(numel(RandMC),1);
sin2T4_bf = zeros(numel(RandMC),1);
for i=RandMC
    progressbar(i/nContours)
    [mnu4Sq{i},sin2T4{i},chi2{i},chi2_ref{i},savefile{i}] = KSN1GridSearch('range',range,...
        'nGridSteps',nGridSteps,...
        'chi2',chi2Str,...
        'DataType',DataType,...
        'freePar','E0 Bkg Norm',...
        'RunList',RunList,...
        'SmartGrid',SmartGrid,...
        'RecomputeFlag','OFF',...
        'RandMC',i,...
        'SysBudget',SysBudget);
    
    chi2tmp   = chi2{i};
    DeltaChi2(i) = chi2tmp(1,1)-min(min(chi2tmp));
    
    % find best fit
    mnu4Sqtmp = mnu4Sq{i};
    sin2T4tmp = sin2T4{i};
    
    [row, col] = find(chi2tmp == min(chi2tmp(:)));
    mnu4Sq_bf(i) =  mnu4Sqtmp(col,row);
    sin2T4_bf(i)  = sin2T4tmp(col,row);
    
end
%%
mnu4Sq_bf = mnu4Sq_bf(RandMC);
sin2T4_bf =sin2T4_bf(RandMC); 

ClosedLog95 =  DeltaChi2>= GetDeltaChi2(0.95,2);
ClosedFrac95 = sum(ClosedLog95)/nContours;
fprintf('%.0f%% C.L. : fraction of significant best fits = %.1f   (%.0f out of %.0f)\n',...
    95,ClosedFrac95*100,sum(ClosedLog95),nContours);
ClosedLog82 =  DeltaChi2>= GetDeltaChi2(0.84,2);
ClosedFrac82 = sum(ClosedLog82)/nContours;
fprintf('%.0f%% C.L. : fraction of significant best fits = %.1f  (%.0f out of %.0f)\n',...
    84,ClosedFrac82*100,sum(ClosedLog82),nContours);
%%
if strcmp(chi2Str,'chi2Stat')
    chi2Label = 'stat. only';
else
    chi2Label = 'stat. and syst.';
end
%% plot only best fits + sensitivity
GetFigure;
%pbf = scatter(sin2T4_bf,mnu4Sq_bf,'filled','MarkerFaceAlpha',0.1,'MarkerFaceColor',rgb('DarkSlateGray'));
%pbf.SizeData=80;
yedge = sort(mnu4Sq_bf);
xedge = sort(sin2T4_bf);
h = histogram2(sin2T4_bf,mnu4Sq_bf,xedge,yedge,'FaceColor','flat','Normalization','probability'); 
PrettyFigureFormat('FontSize',24);
view([0 0 1])
grid off
c = colorbar;
colormap('cool')
c.Label.String = 'Best fit probability';
c.Label.FontSize = get(gca,'FontSize');
set(gca,'YScale','log');
set(gca,'XScale','log');
xlabel('|U_{e4}|^2');
ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
 [mnu4Sq_sensi,sin2T4_sensi,chi2_sensi,chi2_ref_sensi,savefile_sensi] = KSN1GridSearch('range',range,...
        'nGridSteps',nGridSteps,...
        'chi2',chi2Str,...
        'DataType',DataType,...
        'freePar','E0 Bkg Norm',...
        'RunList',RunList,...
        'SmartGrid',SmartGrid,...
        'RecomputeFlag','OFF',...
        'RandMC','OFF');
     PlotArg ={'mnu4Sq',mnu4Sq_sensi,...
        'sin2T4',sin2T4_sensi,...
        'chi2',chi2_sensi,'chi2_ref',chi2_ref_sensi,...
        'CL',CL};

   [pHandle,legStr] =  KSN1ContourPlot(PlotArg{:},'LineStyle','-','PlotSplines','ON','HoldOn','ON','Color','Black');
   pHandle.ZData = 1e3*ones(1e3,1);
xlim([1e-03 0.5]);
ylim([1 94^2]);
leg = legend([h,pHandle],'MC best fits','Sensitivity','EdgeColor',rgb('Silver'),'Location','southwest');
plotdir = strrep(extractBefore(savefile{1},'KSN1'),'results','plots');
MakeDir(plotdir);
plotname = sprintf('%sRandMC_BestFits_%.0f.png',plotdir,range);
print(gcf,plotname,'-dpng','-r450');
%%
close all ;
yedge = sort(mnu4Sq_bf);
xedge = sort(sin2T4_bf);
h = histogram2(sin2T4_bf,mnu4Sq_bf,xedge,yedge,'FaceColor','flat','Normalization','probability'); 
set(gca,'YScale','log');
set(gca,'XScale','log');
view(2)
grid off;
colorbar
%% plot all contours
PlotContoursFlag = 'OFF';
if strcmp(PlotContoursFlag,'ON')
    for i=RandMC
        progressbar(i/nContours)
        titleStr = sprintf('Randomized MC data %.0f (%s) %.0f eV range',i,chi2Label,range);
        SaveAs = sprintf('RandomizedMC/KSN1_GridSearch_KNM1_RandomizedTwin%.0f_%s_%.0feVrange_%s_%.0fnGrid_Grid.png',...
            i,chi2Str,range,strrep(freePar,' ',''),nGridSteps);
        PlotArg ={'mnu4Sq',mnu4Sq{i},...
            'sin2T4',sin2T4{i},...
            'chi2',chi2{i},'chi2_ref',chi2_ref{i},...
            'CL',CL,...
            'titleStr',titleStr,...
            'SaveAs',SaveAs};
        KSN1GridPlot(PlotArg{:},'nInter',1e3);
        close all
    end
end
%%

 %KSN1ContourPlot(PlotArg{:},'LineStyle','-','PlotSplines','ON');
 
 