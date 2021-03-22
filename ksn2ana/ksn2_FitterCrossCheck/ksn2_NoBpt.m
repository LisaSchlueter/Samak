% ksn2 contour on twins without penning trap background slope for fitter
% comparison
% ksn2 calculate chi2 grid search
%% settings that might change
chi2 = 'chi2Stat';
DataType = 'Twin';
nGridSteps = 25;
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
    'FSDFlag','KNM2_0p1eV',...
    'ELossFlag','KatrinT2A20',...
    'AnaFlag','StackPixel',...
    'chi2',chi2,...
    'NonPoissonScaleFactor',NonPoissonScaleFactor,...
    'FSD_Sigma',sqrt(0.0124+0.0025),...
    'TwinBias_FSDSigma',sqrt(0.0124+0.0025),...
    'TwinBias_Q',18573.7,...
    'PullFlag',99,...;%99 = no pull
    'BKG_PtSlope',0,...
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

%%
savedir = sprintf('%sksn2ana/ksn2_FitterCrossCheck/results/',getenv('SamakPath'));
MakeDir(savedir);

S = SterileAnalysis(SterileArg{:});
%S.GridSearch;

S.RunAnaObj.ModelObj.BKG_PtSlope = 3e-06;
S.RunAnaObj.TwinBias_BKG_PtSlope = 3e-06;
S.nGridSteps = 50;
S.LoadGridFile('ExtmNu4Sq','OFF');
S.InterpMode = 'spline';
S.Interp1Grid('RecomputeFlag','ON');
p1 = S.ContourPlot;

Write2Txt('filename',[savedir,'KSN2_contour_Samak_stat_40eV_E0NormBkg'],...
    'Format','dat','variable',[S.sin2T4_contour;S.mNu4Sq_contour],'nCol',2,'variableName','sinT4Sq m4Sq');

S.nGridSteps = 25;
S.RunAnaObj.ModelObj.BKG_PtSlope = 0;
S.RunAnaObj.TwinBias_BKG_PtSlope = 0;
S.LoadGridFile('ExtmNu4Sq','OFF');
S.InterpMode = 'spline';
S.Interp1Grid('RecomputeFlag','ON');
p2 = S.ContourPlot('HoldOn','ON','Color',rgb('Orange'),'LineStyle',':');

S.RunAnaObj.ModelObj.BKG_PtSlope = 0;
S.RunAnaObj.TwinBias_BKG_PtSlope = 3e-06;
S.LoadGridFile('ExtmNu4Sq','OFF');
S.InterpMode = 'spline';
S.Interp1Grid('RecomputeFlag','ON');
p3 = S.ContourPlot('HoldOn','ON','Color',rgb('ForestGreen'),'LineStyle','-.');

Write2Txt('filename',[savedir,'KSN2_contour_Samak_stat_40eV_E0NormBkg_0mucpsPerSModelPtSlope_3mucpsPerSTwinPtSlope'],...
    'Format','dat','variable',[S.sin2T4_contour;S.mNu4Sq_contour],'nCol',2,'variableName','sinT4Sq m4Sq');

xlim([5e-03,0.5])
leg = legend([p1,p2,p3],sprintf('Twin \\alpha = 3 \\mucps/s , Model \\alpha = 3 \\mucps/s '),...
    sprintf('Twin \\alpha = 0 \\mucps/s , Model \\alpha = 0 \\mucps/s '),...
    sprintf('Twin \\alpha = 3 \\mucps/s , Model \\alpha = 0 \\mucps/s '));
leg.Title.String = 'Background time slope from Penning trap';
leg.Title.FontWeight = 'normal';
PrettyLegendFormat(leg);

plotdir = sprintf('%sksn2ana/ksn2_FitterCrossCheck/plots/',getenv('SamakPath'));
MakeDir(plotdir);
plotname = sprintf('%sksn2_BkgPTComparison.png',plotdir);
print(gcf,plotname,'-dpng','-r350');
fprintf('save plot to %s \n',plotname);

