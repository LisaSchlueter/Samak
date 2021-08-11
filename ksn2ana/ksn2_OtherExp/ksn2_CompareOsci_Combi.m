% ksn2 calculate chi2 grid search
%% settings that might change
chi2 = 'chi2CMShape';
DataType = 'Twin';

range = 40;
freePar = 'E0 Norm Bkg';
%% configure RunAnalysis object
if strcmp(chi2,'chi2Stat')
    NonPoissonScaleFactor = 1;
elseif  strcmp(chi2,'chi2CMShape')
    NonPoissonScaleFactor = 1.112;
end
RunAnaArg = {'RunList','KNM2_Prompt',...
    'DataType',DataType,...
    'fixPar',freePar,...%free par
    'SysBudget',40,...
    'fitter','minuit',...
    'minuitOpt','min;migrad',...
    'RadiativeFlag','ON',...
    'FSDFlag','KNM2_0p5eV',...
    'ELossFlag','KatrinT2A20',...
    'AnaFlag','StackPixel',...
    'chi2',chi2,...
    'NonPoissonScaleFactor',NonPoissonScaleFactor,...
    'FSD_Sigma',sqrt(0.0124+0.0025),...
    'TwinBias_FSDSigma',sqrt(0.0124+0.0025),...
    'TwinBias_Q',18573.7,...
    'PullFlag',99,...;%99 = no pull
    'BKG_PtSlope',3*1e-06,...
    'TwinBias_BKG_PtSlope',3*1e-06,...
    'DopplerEffectFlag','FSD'};
A = MultiRunAnalysis(RunAnaArg{:});
%% configure Sterile analysis object
SterileArg = {'RunAnaObj',A,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
    'SmartGrid','OFF',...
    'RecomputeFlag','OFF',...
    'SysEffect','all',...
    'RandMC','OFF',...
    'range',range};

%%
S = SterileAnalysis(SterileArg{:});
if strcmp(DataType,'Real')
    S.nGridSteps = 40;
    S.LoadGridArg = {'ExtmNu4Sq','ON','mNu4SqTestGrid',5};
else
     S.LoadGridArg = {'ExtmNu4Sq','ON','mNu4SqTestGrid',2};
     S.nGridSteps = 30;
end
S.LoadGridFile(S.LoadGridArg{:});
S.Interp1Grid;

FinalSensitivity = 'ON';
[legHandle,legStr,pHandle] = S.ContourPlotOsci('DayaBay','ON','DoubleChooz','ON',...
    'Solar','OFF',...
    'SavePlot','OFF','Style','PRL','BestFit','ON',...
    'FinalSensitivity',FinalSensitivity);

% load combined result for fixed m_nu
savedir_fix = [getenv('SamakPath'),'ksn2ana/ksn2_RunCombination/results/'];
MakeDir(savedir_fix)
savefile_fix = sprintf('%sksn21_Combination_ReAna_%s.mat',savedir_fix,DataType);

% load file
if exist(savefile_fix,'file') 
    dfix = importdata(savefile_fix);
    fprintf('load results from file %s \n',savefile_fix)
else
     fprintf('load results from file %s \n',savefile_fix)
    return
end
mNu4Sq_fix   = dfix.mNu4Sq_contour_12;
sin2T4Sq_fix  = 4*dfix.sin2T4_contour_12.*(1-dfix.sin2T4_contour_12);
mNu4Sq_fix1   = dfix.mNu4Sq_contour_1;
sin2T4Sq_fix1  = 4*dfix.sin2T4_contour_1.*(1-dfix.sin2T4_contour_1);

% include n plot
hold on;
pfix12 = plot(sin2T4Sq_fix,mNu4Sq_fix,':','Color',rgb('Navy'),'LineWidth',2.5);
pfix1 = plot(sin2T4Sq_fix1,mNu4Sq_fix1,'-.','Color',rgb('SkyBlue'),'LineWidth',2.5);

ax = gca;
set(gca,'FontSize',24);
ax.XLabel.FontSize = 26;
ax.YLabel.FontSize = 26;

% legend
legStrCombi = legStr;
if strcmp(FinalSensitivity,'ON')
    legStrCombi{end-2} = [legStr{end-2},' (KSN2)'];
    legStrCombi = {legStrCombi{1:end-3},...
        [legStr{end-2},' (KSN1)'],...
        [legStr{end-2},' (KSN2)'],...
        [legStr{end-2},' (KSN1+2)'],...
        legStrCombi{end-1:end}};
    legHandleAll = [legHandle{1:end-3},...
        pfix1,legHandle{end-2},pfix12,...
        legHandle{end-1},legHandle{end}];
else
    legStrCombi{end} = [legStr{end},' (KSN2)'];
    legStrCombi = {legStrCombi{1:end-1},[legStr{end},' (KSN1)'],legStrCombi{end},[legStr{end},' (KSN1+2)']};
    legHandleAll = [legHandle{1:end-1},pfix1,legHandle{end},pfix12];
end

leg = legend(legHandleAll,legStrCombi);
leg.FontSize=19;
leg.Position(1) = 0.005;
leg.Position(2) = 0.7;
ax.Position(4) = 0.55;
plotname = sprintf('%s_OsciContour_%.2gCL_Combi.pdf',S.DefPlotName,S.ConfLevel);
export_fig(gcf,plotname);
fprintf('save plot to %s \n',plotname)         
print('KSN2_Osci.png','-dpng','-r300');