% investigate impact of different grid types on contour 
nGridSteps = 25;
mNu4SqTestGrid = [1,2,3];
%% configure RunAnalysis object
chi2 = 'chi2Stat';
DataType = 'Twin';
range = 40;
InterpMode = 'spline';
if strcmp(chi2,'chi2Stat')
    NonPoissonScaleFactor = 1;
elseif  strcmp(chi2,'chi2CMShape')
    NonPoissonScaleFactor = 1.112;
end
RunAnaArg = {'RunList','KNM2_Prompt',...
    'DataType',DataType,...
    'fixPar','E0 Norm Bkg',...%free par
    'SysBudget',40,...
    'fitter','minuit',...
    'minuitOpt','min;migrad',...
    'RadiativeFlag','ON',...
    'FSDFlag','KNM2_0p1eV',...
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
    'range',range};
S = SterileAnalysis(SterileArg{:});
%%
for i=1:numel(mNu4SqTestGrid)
S.GridSearch('mNu4SqTestGrid',mNu4SqTestGrid(i),'ExtmNu4Sq','ON');
end

S.GridSearch('mNu4SqTestGrid','OFF','ExtmNu4Sq','ON');
%%
% HoldOn = 'OFF';
% Colors = {'DodgerBlue','Orange','ForestGreen','FireBrick'};
% LineStyles = {'-','-.',':','--'};
% pHandle = cell(numel(nGridSteps),1);
% legStr = cell(numel(nGridSteps),1);
% S.InterpMode = InterpMode;
% for i=1:numel(nGridSteps)
% S.nGridSteps = nGridSteps(i);
% S.LoadGridFile('CheckExtmNu4Sq','ON','CheckLargerN','OFF','CheckSmallerN','OFF');
% S.Interp1Grid('Maxm4Sq',38.2^2);
% pHandle{i} = S.ContourPlot('HoldOn',HoldOn,'Color',rgb(Colors{i}),'LineStyle',LineStyles{i});
% legStr{i} = sprintf('%.0f x %.0f',nGridSteps(i),nGridSteps(i));
% HoldOn = 'ON';
% end
% 
% leg = legend([pHandle{:}],legStr);
% PrettyLegendFormat(leg);
% leg.Title.String = 'Grid size'; leg.Title.FontWeight = 'normal';
% xlim([7e-3,0.5]);
% ylim([1,1.6e3])
% 
% plotdir = [getenv('SamakPath'),'ksn2ana/ksn2_GridSize/plots/'];
% plotname = sprintf('%sksn2_DifferentGrids_%s.png',plotdir,InterpMode);
% MakeDir(plotdir);
% print(plotname,'-dpng','-r300');
% fprintf('save plot to %s \n',plotname);