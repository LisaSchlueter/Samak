% contour with randomized twins
%% settings
CL = 95;
range = 95;%
nGridSteps = 25;
chi2Str = 'chi2CMShape';
DataType = 'Twin';
freePar = 'E0 Bkg Norm';
RunList = 'KNM1';
SmartGrid = 'OFF';
Mode = 'New';
switch Mode
    case 'Old'
        RandMC = [1:151,500:643]*1e3; 
        SysBudget =22;
    case 'New'
        RandMC = [1:155,500:645,901:1101];
        SysBudget =24;
end
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
chi2min_bf   = zeros(numel(RandMC),1);

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
    chi2min_bf(i) =min(chi2tmp(:));
end
%
mnu4Sq_bf = mnu4Sq_bf(RandMC);
sin2T4_bf =sin2T4_bf(RandMC); 
chi2min_bf = chi2min_bf(RandMC);
%% significant best fits

ClosedLog95 =  DeltaChi2>= GetDeltaChi2(0.95,2);
ClosedFrac95 = sum(ClosedLog95)/nContours;
RelErr = @(n,p) sqrt((n*p*(1-p)))/n;

fprintf('%.0f%% C.L. : fraction of significant best fits = %.1f  +- %.1f  (%.0f out of %.0f)\n',...
    95,ClosedFrac95*100,RelErr(nContours,ClosedFrac95)*100,sum(ClosedLog95),nContours);

% ClosedLog82 =  DeltaChi2>= GetDeltaChi2(0.84,2);
% ClosedFrac82 = sum(ClosedLog82)/nContours;
% fprintf('%.0f%% C.L. : fraction of significant best fits = %.1f +- %.1f (%.0f out of %.0f)\n',...
%     84,ClosedFrac82*100,RelErr(nContours,ClosedFrac82)*100,sum(ClosedLog82),nContours);
%% plot only best fits + sensitivity
GetFigure;
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
legend([h,pHandle],'MC best fits',sprintf('Sensitivity %.0f eV range (%.0f%% C.L.)',range,CL),'EdgeColor',rgb('Silver'),'Location','southwest');
plotdir = strrep(extractBefore(savefile{1},'KSN1'),'results','plots');
MakeDir(plotdir);
plotname = sprintf('%sRandMC_BestFits_%.0f_%s.png',plotdir,range,Mode);
print(gcf,plotname,'-dpng','-r450');
%% plot chi2 distribution of best-fits
GetFigure;
hchi2 = histogram(chi2min_bf,'BinWidth',3,...
    'FaceAlpha',1,'FaceColor',rgb('DeepSkyBlue'),'EdgeColor',rgb('SteelBlue'),'Normalization','probability');
hold on;
d = importdata(savefile{1});
dof = d.FitResults_ref.dof;
x = linspace(0,dof*3,1e3);
y = chi2pdf(x,dof);
pchi2 = plot(x,y*hchi2.BinWidth,'Color',rgb('Black'),'LineWidth',2);
PrettyFigureFormat('FontSize',22);
xlabel(sprintf('\\chi^2_{min} (%.0f dof)',dof));
ylabel('Frequency');
leg = legend([hchi2,pchi2],sprintf('%.0f pseudo-experiments',numel(chi2min_bf)),...
                 sprintf('\\chi^2 distribution for %.0f dof',dof),...
                 'EdgeColor',rgb('Silver'));
xlim([0 70]);
plotnameChi2 = sprintf('%sRandMC_BestFitsChi2Dist_%.0f_%s.png',plotdir,range,Mode);
print(gcf,plotnameChi2,'-dpng','-r450');

%% delta chi2 distribution: Kolmogorov-Smirnov Test
close all
GetFigure;
chi2min =sort(chi2min_bf);
Chi2CDFEmp = cumsum(chi2min)./max(cumsum(chi2min));
Chi2CDFTheo = chi2cdf(chi2min,dof);
pEmp =plot(chi2min,Chi2CDFEmp,'-.','LineWidth',2);
hold on;
pTheo = plot(chi2min,Chi2CDFTheo,'-','LineWidth',2);
PrettyFigureFormat('FontSize',22);
xlabel(sprintf('\\chi^2 (%.0f dof)',dof));
ylabel(sprintf('Cumulative probability'));
[h,p,ksstat,cv] = kstest(chi2min,'CDF',[chi2min,Chi2CDFTheo]);
resultsStr = sprintf('KS test: p-value = %.3f',p);
%sprintf('max. difference = %.2f',max(abs(Chi2CDFEmp-Chi2CDFTheo)))
title(resultsStr,'FontWeight','normal','FontSize',get(gca,'FontSize'))
legend([pTheo,pEmp],sprintf('Null hypothesis: \\chi^2 distribution'),...
    sprintf('Empirical distribution'),'EdgeColor',rgb('Silver'),'Location','northwest');
ylim([-0.05 1.05])
plotnameChi2KS = sprintf('%sRandMC_Chi2DistKSTest_%.0f_%s.png',plotdir,range,Mode);
%print(gcf,plotnameChi2KS,'-dpng','-r450');
%% plot all contours - one after the after
if strcmp(chi2Str,'chi2Stat')
    chi2Label = 'stat. only';
else
    chi2Label = 'stat. and syst.';
end
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