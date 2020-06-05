% contour with randomized twins
%% settings
CL = 95;
range = 40;%
nGridSteps = 25;
chi2Str = 'chi2CMShape';
DataType = 'Twin';
freePar = 'E0 Bkg Norm';
RunList = 'KNM1';
SmartGrid = 'OFF';
Mode = 'New_Sterile';
Twin_mNu4Sq   = 1e2;
Twin_sin2T4   = 0.06;

switch Mode
    case 'Old'
        RandMC = [1:151,500:643]*1e3;
        SysBudget =22;
    case 'New'
          SysBudget =24;
        if range== 95      
            RandMC = 1:1500;
        elseif range==65
            RandMC = [1:654,868:931,1219:1500];
        elseif range==40
            RandMC = [1:2000];
        end
    case 'New_Sterile'
        SysBudget =24;    
        RandMC = [1:344,394:456,550:672,728:849];      
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
chi2min_null  = zeros(numel(RandMC),1);

for i=RandMC
    progressbar(i/nContours)
    [mnu4Sq{i},sin2T4{i},chi2{i},chi2_ref{i},savefile{i},FitResults_Null] = KSN1GridSearch('range',range,...
        'nGridSteps',nGridSteps,...
        'chi2',chi2Str,...
        'DataType',DataType,...
        'freePar','E0 Bkg Norm',...
        'RunList',RunList,...
        'SmartGrid',SmartGrid,...
        'RecomputeFlag','OFF',...
        'RandMC',i,...
        'SysBudget',SysBudget,...
        'Twin_mNu4Sq',Twin_mNu4Sq,...
        'Twin_sin2T4',Twin_sin2T4);
    
    chi2tmp   = chi2{i};

    % find best fit
    mnu4Sqtmp = mnu4Sq{i};
    sin2T4tmp = sin2T4{i};
    
    [row, col] = find(chi2tmp == min(chi2tmp(:)));
    mnu4Sq_bf(i) =  mnu4Sqtmp(col,row);
    sin2T4_bf(i)  = sin2T4tmp(col,row);
    chi2min_bf(i) = min(chi2tmp(:));
    chi2min_null(i) = FitResults_Null.chi2min;
    DeltaChi2(i) = chi2min_null(i)-min(min(chi2tmp));
end
%
mnu4Sq_bf = mnu4Sq_bf(RandMC);
sin2T4_bf =sin2T4_bf(RandMC); 
chi2min_bf = chi2min_bf(RandMC);
chi2min_null = chi2min_null(RandMC);
DeltaChi2 = DeltaChi2(RandMC);

d = importdata(savefile{200});
dof = d.FitResults_Null.dof-2;
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
    'RandMC','OFF',...
    'SysBudget',SysBudget);
PlotArg ={'mnu4Sq',mnu4Sq_sensi,...
    'sin2T4',sin2T4_sensi,...
    'chi2',chi2_sensi,'chi2_ref',chi2_ref_sensi,...
    'CL',CL};

[pHandle,legStr] =  KSN1ContourPlot(PlotArg{:},'LineStyle','-','PlotSplines','ON','HoldOn','ON','Color','Black');
pHandle.ZData = 1e3*ones(1e3,1);
xlim([1e-03 0.5]);
if range==95
    ylim([1 94^2]);
elseif range==65
    ylim([1 65^2]);
else
    ylim([1 40^2]);
end
legend([h,pHandle],sprintf('Best fits %.0f pseudo-experiments',numel(RandMC)),sprintf('Sensitivity %.0f eV range (%.0f%% C.L.)',range,CL),'EdgeColor',rgb('Silver'),'Location','southwest');

title(sprintf('Significant best fits (%.2f \\pm %.2f)%% ',...
    ClosedFrac95*100,RelErr(nContours,ClosedFrac95)*100),...
    'FontSize',get(gca,'FontSize'));
plotdir = strrep(extractBefore(savefile{200},'KSN1'),'results','plots');
MakeDir(plotdir);
plotname = sprintf('%sRandMC_BestFits_%.0feV_%.0fsamples.png',plotdir,range,numel(RandMC));
print(gcf,plotname,'-dpng','-r450');
fprintf('save plot to %s \n',plotname);
export_fig(gcf,strrep(plotname,'.png','.pdf'));
%% delta chi2 distribution (best fit chi2): Kolmogorov-Smirnov Test
close all
GetFigure;
chi2min    = sort(chi2min_bf);
Chi2CDFEmp = arrayfun(@(x) sum(chi2min<=x)./numel(chi2min),chi2min);

Chi2CDFTheo = chi2cdf(chi2min,dof);
pEmp =plot(chi2min,Chi2CDFEmp,'-.','LineWidth',2);
hold on;
pTheo = plot(chi2min,Chi2CDFTheo,'-','LineWidth',2);
PrettyFigureFormat('FontSize',22);
xlabel(sprintf('\\chi^2_{min} (%.0f dof)',dof));
ylabel(sprintf('Cumulative probability'));
[h,p,ksstat,cv] = kstest(chi2min,'CDF',[chi2min,Chi2CDFTheo]);
resultsStr = sprintf('KS test: p-value = %.2g',p);
title(resultsStr,'FontWeight','normal','FontSize',get(gca,'FontSize'))
legend([pTheo,pEmp],sprintf(' \\chi^2 distribution with %.0f dof',dof),...
    sprintf(' Empirical distribution (%.0f samples)',numel(chi2min)),...
    'EdgeColor',rgb('Silver'),'Location','southeast');
ylim([-0.05 1.05])
plotnameChi2KS = sprintf('%sRandMC_Chi2DistKSTest_%.0feV_%.0fsamples.png',plotdir,range,numel(RandMC));
print(gcf,plotnameChi2KS,'-dpng','-r450');
fprintf('save plot to %s \n',plotnameChi2KS);
export_fig(gcf,strrep(plotnameChi2KS,'.png','.pdf'));
%% plot chi2 distribution of best-fits
GetFigure;
hchi2 = histogram(chi2min_bf,'BinWidth',3,...
    'FaceAlpha',1,'FaceColor',rgb('DeepSkyBlue'),'EdgeColor',rgb('SteelBlue'),'Normalization','probability');
hold on;
x = linspace(0,dof*3,1e3);
y = chi2pdf(x,dof);
pchi2 = plot(x,y*hchi2.BinWidth,'Color',rgb('Black'),'LineWidth',2);
PrettyFigureFormat('FontSize',22);
xlabel(sprintf('\\chi^2_{min} (%.0f dof)',dof));
ylabel('Frequency');
title(resultsStr,'FontWeight','normal','FontSize',get(gca,'FontSize'))
leg = legend([hchi2,pchi2],sprintf('%.0f pseudo-experiments',numel(chi2min_bf)),...
                 sprintf('\\chi^2 distribution for %.0f dof',dof),...
                 'EdgeColor',rgb('Silver'));
xlim([0 70]);
plotnameChi2 = sprintf('%sRandMC_BestFitsChi2Dist_%.0feV_%.0fsamples.png',plotdir,range,numel(RandMC));
print(gcf,plotnameChi2,'-dpng','-r450');
fprintf('save plot to %s \n',plotnameChi2);
export_fig(gcf,strrep(plotnameChi2,'.png','.pdf'));

%% DeltaChi2 distribution (best fit - null chi2)
GetFigure;
PlotDeltaChi2 = sort(DeltaChi2);
DeltaChi2CDF = arrayfun(@(x) sum(PlotDeltaChi2<=x)./numel(PlotDeltaChi2),PlotDeltaChi2);
% calculate 95 quantile: interpolation
[DeltaChi2CDFquantile,ia] = unique(DeltaChi2CDF);
DeltaChi2CrApprox = interp1(DeltaChi2CDFquantile,PlotDeltaChi2(ia),0.95,'lin');%quantile(PlotDeltaChi2,0.95);% PlotDeltaChi2(find(abs(DeltaChi2CDF-0.95)==min(abs(DeltaChi2CDF-0.95)),1));

% plot
DeltaChi2CDFTheo = chi2cdf(PlotDeltaChi2,2);
xInter = linspace(0,10,1e2);
DeltaChi2CrTheo = interp1(chi2cdf(xInter,2),xInter,0.95,'spline');
p95 = plot(linspace(0,20,10),0.95*ones(10,1),'-','LineWidth',1.5,'Color',rgb('Silver'));
hold on;
pchi2Theo = plot(PlotDeltaChi2,DeltaChi2CDFTheo,'-','LineWidth',2.5,'Color',rgb('Orange'));
pchi2 = plot(PlotDeltaChi2,DeltaChi2CDF,'-.','LineWidth',2.5,'Color',rgb('DodgerBlue'));

PrettyFigureFormat('FontSize',22);
xlabel(sprintf('\\Delta \\chi^2'));
ylabel(sprintf('Cumulative probability'));
xlim([0 10]);
ylim([-0.02 1.02]);
legend([p95,pchi2Theo,pchi2,],sprintf('95%% quantile'),...
    sprintf('\\chi^2 distribution for 2 dof          \\Delta\\chi^2_{crit.} = %.2f',DeltaChi2CrTheo),...
    sprintf('Empirical cdf (%.0f samples) \\Delta\\chi^2_{crit.} = %.2f',numel(PlotDeltaChi2),DeltaChi2CrApprox),...
    'EdgeColor',rgb('Silver'),'Location','southeast');
t = title(sprintf('\\Delta\\chi^2 = \\chi^2_{null} - \\chi^2_{min} '),'FontWeight','normal','FontSize',get(gca,'FontSize'));
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
plotnameChi2crit = sprintf('%sRandMC_DeltaChi2Crit_%.0feV_%.0fsamples.png',plotdir,range,numel(RandMC));
print(gcf,plotnameChi2crit,'-dpng','-r450');
fprintf('save plot to %s \n',plotnameChi2crit);
export_fig(gcf,strrep(plotnameChi2crit,'.png','.pdf'));
%% calculate quantile of you best fit (delta chi2 = 1.58)
DeltaChi2_bf = 1.58;
NumLargerDeltaChi2 = sum(PlotDeltaChi2>=DeltaChi2_bf);
FracLargerDeltaChi2 = NumLargerDeltaChi2/numel(PlotDeltaChi2);

fprintf('Empirical: %.0f out of %.0f have Delta chi2  larger or equal compared to best fit (%.2f) -> p-value = %.2f \n',...
    NumLargerDeltaChi2,numel(PlotDeltaChi2),DeltaChi2_bf,FracLargerDeltaChi2)
fprintf('Empirical: p-value from interpolation instead of fraction: %.2f \n',1-interp1(PlotDeltaChi2(ia),DeltaChi2CDFquantile,DeltaChi2_bf,'lin'));
fprintf('p-value from expected chi2cdf (2 dof) = %.3f \n',1-chi2cdf(DeltaChi2_bf,2));
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