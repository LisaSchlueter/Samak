% Test SterileAnalysis class
% Lisa, May2020
% plot with neutrino mass as
% nuissance parameter free
% nuissance parameter + pull 12 (Mainz & Troitsk)
% nuissance parameter + pull 13 (E0 1 eV)
% nuissance parameter + pull 14 (E0 2 eV)
%% settings for runanalysis
DataType = 'Real';
%%
RunAnaArg = {'RunList','KNM1',...
    'fixPar','mNu E0 Norm Bkg',...
    'DataType',DataType,...
    'FSDFlag','SibilleFull',...
    'ELossFlag','KatrinT2',...
    'AnaFlag','StackPixel',...
    'chi2','chi2Stat',...
    'ROIFlag','Default',...
    'SynchrotronFlag','ON',...
    'AngularTFFlag','OFF',...
    'ISCSFlag','Edep',...
    'TwinBias_Q',18573.73,...
    'SysBudget',24,...
    'pullFlag',99,...
    'NonPoissonScaleFactor',1};

T = MultiRunAnalysis(RunAnaArg{:});
T.chi2 = 'chi2CMShape';
%% settings sterile class
SterileArg = {'RunAnaObj',T,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
    'nGridSteps',50,...
    'SmartGrid','OFF',...
    'RecomputeFlag','OFF',...
    'SysEffect','all',...
    'RandMC','OFF',...
    'range',40};

S = SterileAnalysis(SterileArg{:});

%% 2 load a chi2 map 1
S.InterpMode = 'spline'; % waring: if contour is closed, spline interp sometimes sensitive to artefacts! Switch to "lin" in this case
S.RunAnaObj.pullFlag = 99;
S.RunAnaObj.fixPar = 'E0 Norm Bkg'; S.RunAnaObj.InitFitPar;
S.LoadGridFile('CheckSmallerN','ON','CheckLargerN','ON'); % if CheckSmallerN also look for grid with more/less nGridSteps
S.Interp1Grid('RecomputeFlag','ON');% interpolate chi2 map -> nicer appearance of all plots. some
pfix =S.ContourPlot('HoldOn','OFF','Color',rgb('LimeGreen'),'LineStyle','-');

S.InterpMode = 'lin';
S.RunAnaObj.fixPar = 'mNu E0 Norm Bkg'; S.RunAnaObj.InitFitPar;
S.RunAnaObj.pullFlag = 99;
S.LoadGridFile('CheckSmallerN','ON','CheckLargerN','ON'); % if CheckSmallerN also look for grid with more/less nGridSteps
S.Interp1Grid('RecomputeFlag','ON');% interpolate chi2 map -> nicer appearance of all plots. some
pfree = S.ContourPlot('HoldOn','ON','Color',rgb('DodgerBlue'),'LineStyle','-.');

S.InterpMode = 'spline';
S.RunAnaObj.pullFlag = 12;
S.LoadGridFile('CheckSmallerN','ON','CheckLargerN','ON'); % if CheckSmallerN also look for grid with more/less nGridSteps
S.Interp1Grid('RecomputeFlag','ON');% interpolate chi2 map -> nicer appearance of all plots. some
ppull1 = S.ContourPlot('HoldOn','ON','Color',rgb('Navy'),'LineStyle',':');

S.InterpMode = 'lin';
S.RunAnaObj.pullFlag = 13;
S.LoadGridFile('CheckSmallerN','ON','CheckLargerN','ON'); % if CheckSmallerN also look for grid with more/less nGridSteps
S.Interp1Grid('RecomputeFlag','ON');% interpolate chi2 map -> nicer appearance of all plots. some
ppull2 = S.ContourPlot('HoldOn','ON','Color',rgb('Orange'),'LineStyle','--');

leg = legend(...
    sprintf('KATRIN KSN1 %.0f%% C.L. - {\\itm}_\\nu^2 = 0 eV^2',S.ConfLevel),...
    sprintf('KATRIN KSN1 %.0f%% C.L. - {\\itm}_\\nu^2 free',S.ConfLevel),...
    sprintf('KATRIN KSN1 %.0f%% C.L. - {\\itm}_\\nu^2 free - \\sigma({\\itm}_\\nu^2) = 1.94 eV^2',S.ConfLevel),...
    sprintf('KATRIN KSN1 %.0f%% C.L. - {\\itm}_\\nu^2 free - \\sigma({\\itE}_0) = 1 eV',S.ConfLevel));
legend boxoff;

% style
PRLFormat;
set(gca,'FontSize',20);
t = get(gca,'title');
t.FontSize = 20;
xlim([5e-03 0.4])
ylim([1 2e3])

%% save
plotname = strrep(S.DefPlotName,'_pull13','ComparePulls.png');
print(gcf,plotname,'-dpng','-r300');
