% ksn2 calculate chi2 grid search
%% settings that might change
chi2 = 'chi2CMShape';
DataType = 'Real';
nGridSteps = 40;
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
    'nGridSteps',nGridSteps,...
    'SmartGrid','OFF',...
    'RecomputeFlag','OFF',...
    'SysEffect','all',...
    'RandMC','OFF',...
    'range',range,...
    'LoadGridArg',{'ExtmNu4Sq','ON','mNu4SqTestGrid',5}};

%%
S = SterileAnalysis(SterileArg{:});
S.LoadGridFile(S.LoadGridArg{:});
S.Interp1Grid;

%%
 [legHandle,legStr] = S.ContourPlotOsci('DayaBay','ON','DoubleChooz','ON','SavePlot','ON','Style','PRL','BestFit','ON');

% %% load combi- free mNuSq
% savedir = [getenv('SamakPath'),'SterileAnalysis/GridSearchFiles/Combi/',DataType,'/'];
% MakeDir(savedir)
% savename_free = sprintf('%sKSN12Combi_ReAna_GridSearch_%s_%s_Uniform_%s_%.0fnGrid.mat',...
%     savedir,DataType,'mNuE0BkgNorm','chi2CMShape',30);
% 
% % load file
% if exist(savename_free,'file') 
%     dfree = importdata(savename_free);
%     fprintf('load results from file %s \n',savename_free)
% else
%      fprintf('load results from file %s \n',savename_free)
%     return
% end
% 
% mNu4Sq_free    = dfree.mNu4Sq_contour;
% sin2T4Sq_free  = 4*dfree.sin2T4_contour.*(1-dfree.sin2T4_contour);

% load combi fixed m_nu
savedir_fix = [getenv('SamakPath'),'ksn2ana/ksn2_RunCombination/results/'];
MakeDir(savedir)
savefile_fix = sprintf('%sksn21_Combination_ReAna.mat',savedir_fix);

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

%% plot
hold on;
pfix = plot(sin2T4Sq_fix,mNu4Sq_fix,':','Color',rgb('Navy'),'LineWidth',2.5);

%%
%pfree = plot(sin2T4Sq_free,mNu4Sq_free,':','Color',rgb('DarkRed'),'LineWidth',2.5);

legStrCombi = legStr;
legStrCombi{end} = [legStr{end},' (KNM-2)'];
legStrCombi = {legStrCombi{:},[legStr{end},' (KNM-1&2)']};
%%
leg = legend([legHandle{:},pfix],legStrCombi);
 plotname = sprintf('%s_OsciContour_%.2gCL_Combi.pdf',S.DefPlotName,S.ConfLevel);
export_fig(gcf,plotname);
                  