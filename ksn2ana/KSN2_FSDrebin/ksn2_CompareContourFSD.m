% compare contour with different FSDs 
% twin, stat, 40 eV
% ksn2 contour on twins without penning trap background slope for fitter
% comparison
% ksn2 calculate chi2 grid search
%% settings that might change
chi2 = 'chi2Stat';
DataType = 'Twin';
nGridSteps = 25;
range = 40;
FSDFlag = 'KNM2_0p1eV';
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
    'FSDFlag',FSDFlag,...
    'ELossFlag','KatrinT2A20',...
    'AnaFlag','StackPixel',...
    'chi2',chi2,...
    'NonPoissonScaleFactor',NonPoissonScaleFactor,...
    'FSD_Sigma',sqrt(0.0124+0.0025),...
    'TwinBias_FSDSigma',sqrt(0.0124+0.0025),...
    'TwinBias_Q',18573.7,...
    'PullFlag',99,...;%99 = no pull
    'BKG_PtSlope',3e-06,...
    'TwinBias_BKG_PtSlope',3e-06,...
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
S.GridSearch;

%% get number of bin
A.ModelObj.TTFSD = 'KNM2';
A.ModelObj.LoadFSD;
nBin0 = numel(A.ModelObj.TTexE);

A.ModelObj.TTFSD = 'KNM2_0p1eV';
A.ModelObj.LoadFSD;
nBin1 = numel(A.ModelObj.TTexE);

A.ModelObj.TTFSD = 'KNM2_0p1eV_cut40eV';
A.ModelObj.LoadFSD;
nBin2 = numel(A.ModelObj.TTexE);

A.ModelObj.TTFSD = 'KNM2_0p1eV_cut50eV';
A.ModelObj.LoadFSD;
nBin3 = numel(A.ModelObj.TTexE);

A.ModelObj.TTFSD = 'KNM2_0p5eV';
A.ModelObj.LoadFSD;
nBin4 = numel(A.ModelObj.TTexE);

A.ModelObj.TTFSD = 'BlindingKNM2';
A.ModelObj.LoadFSD;
nBin5 = numel(A.ModelObj.TTexE);
%% plot

A.FSDFlag = 'KNM2';
S.LoadGridFile('CheckLargerN','OFF'); 
dtmp = importdata(S.GridFilename); tCpuHour0 = dtmp.tCpuHour;
S.Interp1Grid('maxM4Sq',36^2)
p0 = S.ContourPlot('HoldOn','OFF','Color',rgb('DodgerBlue'),'LineStyle','-');

A.FSDFlag = 'KNM2_0p1eV';
S.LoadGridFile('CheckLargerN','OFF');
dtmp = importdata(S.GridFilename); tCpuHour1 = dtmp.tCpuHour;
S.Interp1Grid('maxM4Sq',36^2)
p1 = S.ContourPlot('HoldOn','ON','Color',rgb('Orange'),'LineStyle','-.');

% A.FSDFlag = 'KNM2_0p1eV_cut40eV';
% S.LoadGridFile('CheckLargerN','OFF');
%dtmp = importdata(S.GridFilename); tCpuHour2 = dtmp.tCpuHour;
% S.Interp1Grid('maxM4Sq',36^2)
% p2 = S.ContourPlot('HoldOn','ON','Color',rgb('ForestGreen'),'LineStyle',':');

% A.FSDFlag = 'KNM2_0p1eV_cut50eV';
% S.LoadGridFile('CheckLargerN','OFF');
%dtmp = importdata(S.GridFilename); tCpuHour3 = dtmp.tCpuHour;
% S.Interp1Grid('maxM4Sq',36^2)
% p3 = S.ContourPlot('HoldOn','ON','Color',rgb('FireBrick'),'LineStyle','--');

A.FSDFlag = 'KNM2_0p5eV';
S.LoadGridFile('CheckLargerN','OFF');
dtmp = importdata(S.GridFilename); tCpuHour4 = dtmp.tCpuHour;
S.Interp1Grid('maxM4Sq',36^2)
p4 = S.ContourPlot('HoldOn','ON','Color',rgb('LimeGreen'),'LineStyle',':');

A.ModelObj.TTFSD = 'BlindingKNM2';
S.LoadGridFile('CheckLargerN','OFF');
dtmp = importdata(S.GridFilename); tCpuHour4 = dtmp.tCpuHour;
S.Interp1Grid('maxM4Sq',36^2)
p5 = S.ContourPlot('HoldOn','ON','Color',rgb('FireBrick'),'LineStyle','--');

%%
leg = legend([p0,p1,p4,p5],...
    sprintf('KNM2 default (%.0f bins)',nBin0),...
    sprintf('KNM2 0.1eV rebin (%.0f bins)',nBin1),...
    sprintf('KNM2 0.5eV rebin (%.0f bins)',nBin4),...
    sprintf('Blinding KNM2 (%.0f bins)',nBin5),...
    'Location','southwest');
%   sprintf('KNM2 0.1eV rebin - cut 40 eV (%.0f bins)',nBin2),...
%   sprintf('KNM2 0.1eV rebin - cut 50 eV (%.0f bins)',nBin3),...
PrettyLegendFormat(leg);
xlim([6e-03,0.5])
%% save
plotdir = sprintf('%sksn2ana/ksn2_FSDrebin/plots/',getenv('SamakPath'));
MakeDir(plotdir);
plotname = sprintf('%sksn2_CompareContourFSD.png',plotdir);
print(gcf,plotname,'-dpng','-r350');
fprintf('save plot to %s \n',plotname);

