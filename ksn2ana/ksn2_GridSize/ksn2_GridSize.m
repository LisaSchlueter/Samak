% investigate impact of grid size on contour
freePar = 'mNu E0 Norm Bkg';
if contains(freePar,'mNu')
    nGridSteps = [30,40];
    chi2 = 'chi2CMShape';
     LoadGridArg = {'CheckLargerN','OFF','CheckSmallerN','OFF','IgnoreKnm2FSDbinning','ON','mNu4SqTestGrid',5};
else
    nGridSteps = [50,30,25,10];
    chi2 = 'chi2Stat';
    LoadGridArg = {'CheckLargerN','OFF','CheckSmallerN','OFF','IgnoreKnm2FSDbinning','ON'};
end
%% configure RunAnalysis object

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
    'fixPar',freePar,...%free par
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
    'nGridSteps',nGridSteps(1),...
    'SmartGrid','OFF',...
    'RecomputeFlag','OFF',...
    'SysEffect','all',...
    'RandMC','OFF',...
    'range',range};
S = SterileAnalysis(SterileArg{:});
%%
HoldOn = 'OFF';
Colors = {'DodgerBlue','Orange','LimeGreen','FireBrick'};
LineStyles = {'-',':','--','-.'};
pHandle = cell(numel(nGridSteps),1);
legStr = cell(numel(nGridSteps),1);
S.InterpMode = InterpMode;
for i=1:numel(nGridSteps)
    if nGridSteps(i)==10 || contains(freePar,'mNu')
        ExtmNu4Sq = 'OFF';
    else
        ExtmNu4Sq = 'ON';
    end
S.nGridSteps = nGridSteps(i);
S.LoadGridFile(LoadGridArg{:},'ExtmNu4Sq',ExtmNu4Sq);
S.Interp1Grid('Maxm4Sq',38.2^2);
pHandle{i} = S.ContourPlot('HoldOn',HoldOn,'Color',rgb(Colors{i}),'LineStyle',LineStyles{i});
legStr{i} = sprintf('%.0f x %.0f',nGridSteps(i),nGridSteps(i));
HoldOn = 'ON';
end

leg = legend([pHandle{:}],legStr);
PrettyLegendFormat(leg);
leg.Title.String = 'Grid size'; leg.Title.FontWeight = 'normal';
xlim([7e-3,0.5]);
ylim([1,1.6e3])
%%
plotdir = [getenv('SamakPath'),'ksn2ana/ksn2_GridSize/plots/'];
plotname = sprintf('%sksn2_GridSize_%s_%s.png',plotdir,InterpMode,strrep(freePar,' ',''));
MakeDir(plotdir);
print(plotname,'-dpng','-r300');
fprintf('save plot to %s \n',plotname);
