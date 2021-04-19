% compare contour with different FSDs 
% real, stat & syst, 40 eV
%% settings that might change
chi2 = 'chi2CMShape';
DataType = 'Real';
nGridSteps = 30;
range = 40;
FSDFlags = {'KNM2','KNM2_0p1eV','KNM2_0p5eV'};
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
    'FSDFlag',FSDFlags{1},...
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
    'range',range,...
    'LoadGridArg',{'mNu4SqTestGrid',5,'ExtmNu4Sq','ON'}};
S = SterileAnalysis(SterileArg{:});

%% get number of bins
nBins = zeros(numel(FSDFlags),1);
for i=1:numel(FSDFlags)
    
    A.ModelObj.TTFSD = FSDFlags{i};
    A.ModelObj.LoadFSD;
    nBins(i) = numel(A.ModelObj.TTexE);
    
end
%% plot
S.InterpMode = 'spline';
pHandles = cell(numel(FSDFlags),1);
legStr = cell(numel(FSDFlags),1);
HoldOn = 'ON';
Colors = {'DodgerBlue','Orange','FireBrick','ForestGreen'};
LineStyle = {'-','-.','--',':','-'};
HoldOn = 'OFF';
tCpuHour = zeros(numel(FSDFlags),1);

for i=1:numel(FSDFlags)
A.FSDFlag = FSDFlags{i};
S.LoadGridFile(S.LoadGridArg{:}); 
dtmp = importdata(S.GridFilename(S.LoadGridArg{:})); tCpuHour(i) = dtmp.tCpuHour;
S.Interp1Grid('maxM4Sq',40^2)
pHandles{i} = S.ContourPlot('HoldOn',HoldOn,'Color',rgb(Colors{i}),'LineStyle',LineStyle{i});
legStr{i} = sprintf('%s (%.0f bins)',strrep(FSDFlags{i},'_',' '),nBins(i));
HoldOn = 'ON';
end


%%
leg = legend([pHandles{:}],legStr,'Location','southwest');
PrettyLegendFormat(leg);
xlim([4e-03,0.5])
ylim([1 1600]);

%% save
plotdir = sprintf('%sksn2ana/ksn2_FSDrebin/plots/',getenv('SamakPath'));
MakeDir(plotdir);
plotname = sprintf('%sksn2_CompareContourFSD_Data.png',plotdir);
print(gcf,plotname,'-dpng','-r350');
fprintf('save plot to %s \n',plotname);

