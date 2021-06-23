% ksn2 calculate chi2 grid search
%% settings that might change
chi2 = 'chi2CMShape';
DataType = 'Twin';
nGridSteps = 25;
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

%%
S = SterileAnalysis(SterileArg{:});
if strcmp(S.RunAnaObj.DataType,'Twin')
S.InterpMode = 'spline';
else
    S.InterpMode = 'lin';
end

%% load and plot ksn1
S.nGridSteps = 50;
S.RunAnaObj.DataSet = 'Knm1';
S.RunAnaObj.RunData.RunName = 'KNM1';
S.RunAnaObj.ELossFlag  = 'KatrinT2';
S.RunAnaObj.AngularTFFlag ='OFF';
S.RunAnaObj.SysBudget = 24;

% stat
S.RunAnaObj.chi2 = 'chi2Stat';
S.LoadGridFile('CheckLarger','ON');
S.Interp1Grid('RecomputeFlag','ON');
p1stat = S.ContourPlot('BestFit','OFF','CL',95,'HoldOn','OFF','Color',rgb('Orange'),'LineStyle','-.');

% stat and syst
S.RunAnaObj.chi2 = 'chi2CMShape';
S.LoadGridFile('CheckLarger','ON');
S.Interp1Grid('RecomputeFlag','ON');
p1tot = S.ContourPlot('BestFit','OFF','CL',95,'HoldOn','ON','Color',rgb('FireBrick'),'LineStyle','-');


% load ans plot ksn2 
S.nGridSteps = 25;
S.RunAnaObj.DataSet = 'Knm2';
S.RunAnaObj.RunData.RunName = 'KNM2_Prompt';
S.RunAnaObj.ELossFlag  = 'KatrinT2A20';
S.RunAnaObj.AngularTFFlag ='ON';
S.RunAnaObj.SysBudget = 40;

% stat
S.RunAnaObj.chi2 = 'chi2Stat';
S.LoadGridFile('CheckLarger','ON');
S.Interp1Grid('RecomputeFlag','ON');
[p2Stat,sinMin] = S.ContourPlot('BestFit','OFF','CL',95,'HoldOn','ON','Color',rgb('SkyBlue'),'LineStyle','-.');

% stat and syst
S.RunAnaObj.chi2 = 'chi2CMShape';
S.LoadGridFile('CheckLarger','ON');
S.Interp1Grid('RecomputeFlag','ON');
[p2tot,sinMin] = S.ContourPlot('BestFit','OFF','CL',95,'HoldOn','ON','Color',rgb('DodgerBlue'),'LineStyle','-');

xlim([sinMin*0.5,0.5]);

%% legend
leg = legend([p1stat,p1tot,p2Stat,p2tot],'KSN-1 (stat. only)','KSN-1 (stat. and syst.)','KSN-2 (stat. only)','KSN-2 (stat. and syst.)');
PrettyLegendFormat(leg);
if ~contains(freePar,'mNu')
t = title(sprintf('%s , {\\itm}_\\nu^2 = 0 eV^2',extractBefore(S.GetPlotTitle,'(')),'FontWeight','normal','FontSize',get(gca,'FontSize')-2);
end

%% save plot
plotname = sprintf('%sCompareKSN12_%s_%.0feV_%s.png',S.DefPlotName,A.DataType,S.range,A.chi2);
print(gcf,plotname,'-r300','-dpng');
fprintf('save plot to %s \n',plotname);