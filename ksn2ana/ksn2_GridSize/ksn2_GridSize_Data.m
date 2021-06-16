% investigate impact of grid size on contour 
nGridSteps = [30,40,50];

%% configure RunAnalysis object
chi2 = 'chi2CMShape';
DataType = 'Real';
range = 40;
InterpMode = 'spline';
freePar = 'E0 Norm Bkg';
PullFlag = 99;%26;

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
    'PullFlag',PullFlag,...;%99 = no pull
    'BKG_PtSlope',3*1e-06,...
    'TwinBias_BKG_PtSlope',3*1e-06,...
    'DopplerEffectFlag','FSD'};
A = MultiRunAnalysis(RunAnaArg{:});
%% configure Sterile analysis object
SterileArg = {'RunAnaObj',A,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
    'nGridSteps',nGridSteps(1),...
    'SmartGrid','OFF',...
    'RecomputeFlag','OFF',...
    'SysEffect','all',...
    'RandMC','OFF',...
    'range',range};
S = SterileAnalysis(SterileArg{:});
%%
HoldOn = 'OFF';
Colors = {'DodgerBlue','Orange','FireBrick','FireBrick'};
MarkerStyle = {'x','*','o'};

LineStyles = {'-',':','--','-.'};
pHandle = cell(numel(nGridSteps),1);
legStr = cell(numel(nGridSteps),1);
S.InterpMode = InterpMode;
for i=1:numel(nGridSteps)
    if nGridSteps(i)==30 || ~contains(freePar,'mNu')
        ExtmNu4Sq = 'ON';
    else
        ExtmNu4Sq = 'OFF';
    end
S.nGridSteps = nGridSteps(i);
S.LoadGridFile('CheckLargerN','OFF','CheckSmallerN','OFF','ExtmNu4Sq',ExtmNu4Sq,'mNu4SqTestGrid',5);
S.Interp1Grid('Maxm4Sq',40^2);
pHandle{i} = S.ContourPlot('BestFit','ON','HoldOn',HoldOn,'Color',rgb(Colors{i}),'LineStyle',LineStyles{i},'MarkerStyle',MarkerStyle{i});
legStr{i} = sprintf('%.0f x %.0f',nGridSteps(i),nGridSteps(i));
HoldOn = 'ON';
end

leg = legend([pHandle{:}],legStr);
PrettyLegendFormat(leg);
leg.Title.String = 'Grid size'; leg.Title.FontWeight = 'normal';
xlim([3e-3,0.5]);
ylim([1,1.6e3])
%%
plotdir = [getenv('SamakPath'),'ksn2ana/ksn2_GridSize/plots/'];
plotname = sprintf('%sksn2_GridSize_Data_%s_%s_%s.png',plotdir,strrep(freePar,' ',''),chi2,InterpMode);
if PullFlag~=99
  plotname = strrep(plotname,'.png',sprintf('_pull%.0f.png',PullFlag));
end
MakeDir(plotdir);
print(plotname,'-dpng','-r300');
fprintf('save plot to %s \n',plotname);
