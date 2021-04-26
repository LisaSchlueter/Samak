% ksn2 chi2 grid - non physical parameter space 4 quadrants
% ksn2 plot chi2 grid search
% 4 quadrants
%% settings that might change
chi2 = 'chi2CMShape';
DataType = 'Real';
nGridSteps = 30;
range = 40;

%% configure RunAnalysis object
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
SterileArg = {... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
    'nGridSteps',nGridSteps,...
    'SmartGrid','OFF',...
    'RecomputeFlag','OFF',...
    'SysEffect','all',...
    'RandMC','OFF',...
    'range',range,...
    'LoadGridArg',{'ExtmNu4Sq','ON','mNu4SqTestGrid',5},...
    'InterpMode','spline'};
S = SterileAnalysis('RunAnaObj',A,SterileArg{:});
%%
S.LoadGridFile(S.LoadGridArg{:});
S.Interp1Grid;
p1 = S.ContourPlot;

%
Chi2Crit = 6.19; % +- 0.12 (from Martin Slezak, H0, 5000 contours)
myCL = GetCL(Chi2Crit);
p2 = S.ContourPlot('CL',myCL,'HoldOn','ON','Color',rgb('Skyblue'),'LineStyle','-.');

xlim([5e-03, 0.5]);
ylim([2, 1600]);

leg = legend([p1,p2],sprintf('\\chi^2_{crit.} = 5.99 (Wilk`s theorem)'),sprintf('\\chi^2_{crit.} = %.2f',Chi2Crit));
PrettyLegendFormat(leg);

savedir = [getenv('SamakPath'),'ksn2ana/ksn2_WilksTheorem/plots/'];
MakeDir(savedir);
pltname = sprintf('%sksn2_WT_ContourShift_Chi2Crit%.2f.pdf',savedir,Chi2Crit);
export_fig(pltname);
fprintf('save olot to %s \n',pltname)
%%
function cl = GetCL(chi2crit)
x = linspace(80,99.99,1e3);
chi2 = GetDeltaChi2(x,2);
cl = interp1(chi2,x,chi2crit);
end
