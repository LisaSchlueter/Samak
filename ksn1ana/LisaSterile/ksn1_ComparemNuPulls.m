% Test SterileAnalysis class
% Lisa, May2020
% plot with neutrino mass as
% nuissance parameter free
% nuissance parameter + 1 eV^2 mnu pull
% nuissance parameter + 2 eV^2 mnu pull
% nuissance parameter + 3 eV^2 mnu pull
%% settings for runanalysis
DataType = 'Twin';
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

 switch S.RunAnaObj.DataType
     case 'Real'
         pullInt = 'lin';
     case 'Twin'
         pullInt = 'spline';
 end

%% 2 load a chi2 map 1
S.InterpMode = 'spline'; % waring: if contour is closed, spline interp sometimes sensitive to artefacts! Switch to "lin" in this case
S.RunAnaObj.pullFlag = 99;
S.RunAnaObj.fixPar = 'E0 Norm Bkg'; S.RunAnaObj.InitFitPar;
S.LoadGridFile('CheckSmallerN','ON','CheckLargerN','ON'); % if CheckSmallerN also look for grid with more/less nGridSteps
S.Interp1Grid('RecomputeFlag','ON');% interpolate chi2 map -> nicer appearance of all plots. some
pfix =S.ContourPlot('HoldOn','OFF','Color',rgb('Black'),'LineStyle','-');
%pfix.LineWidth = 3;

S.InterpMode = pullInt;
S.RunAnaObj.fixPar = 'mNu E0 Norm Bkg'; S.RunAnaObj.InitFitPar;
S.RunAnaObj.pullFlag = 99;
S.LoadGridFile('CheckSmallerN','ON','CheckLargerN','ON'); % if CheckSmallerN also look for grid with more/less nGridSteps
S.Interp1Grid('RecomputeFlag','ON');% interpolate chi2 map -> nicer appearance of all plots. some
pfree = S.ContourPlot('HoldOn','ON','Color',rgb('SlateGray'),'LineStyle','-');
%pfree.LineWidth = 3;

S.InterpMode = pullInt;
S.RunAnaObj.pullFlag = 18;
S.LoadGridFile('CheckSmallerN','ON','CheckLargerN','ON'); % if CheckSmallerN also look for grid with more/less nGridSteps
S.Interp1Grid('RecomputeFlag','ON');% interpolate chi2 map -> nicer appearance of all plots. some
p0p5eV = S.ContourPlot('HoldOn','ON','Color',rgb('IndianRed'),'LineStyle',':');

S.InterpMode =pullInt;
S.RunAnaObj.pullFlag = 15;
S.LoadGridFile('CheckSmallerN','ON','CheckLargerN','ON'); % if CheckSmallerN also look for grid with more/less nGridSteps
S.Interp1Grid('RecomputeFlag','ON');% interpolate chi2 map -> nicer appearance of all plots. some
p1eV = S.ContourPlot('HoldOn','ON','Color',rgb('Orange'),'LineStyle','--');

S.InterpMode = pullInt;
S.RunAnaObj.pullFlag = 16;
S.LoadGridFile('CheckSmallerN','ON','CheckLargerN','ON'); % if CheckSmallerN also look for grid with more/less nGridSteps
S.Interp1Grid('RecomputeFlag','ON');% interpolate chi2 map -> nicer appearance of all plots. some
p2eV = S.ContourPlot('HoldOn','ON','Color',rgb('DodgerBlue'),'LineStyle','-.');

S.InterpMode = pullInt;
S.RunAnaObj.pullFlag = 17;
S.LoadGridFile('CheckSmallerN','ON','CheckLargerN','ON'); % if CheckSmallerN also look for grid with more/less nGridSteps
S.Interp1Grid('RecomputeFlag','ON');% interpolate chi2 map -> nicer appearance of all plots. some
p3eV = S.ContourPlot('HoldOn','ON','Color',rgb('LimeGreen'),'LineStyle',':');

% S.InterpMode = 'lin';
% S.RunAnaObj.pullFlag = 13;
% S.LoadGridFile('CheckSmallerN','ON','CheckLargerN','ON'); % if CheckSmallerN also look for grid with more/less nGridSteps
% S.Interp1Grid('RecomputeFlag','ON');% interpolate chi2 map -> nicer appearance of all plots. some
% pE0 = S.ContourPlot('HoldOn','ON','Color',rgb('Silver'),'LineStyle','--');

leg = legend([pfix,p0p5eV,p1eV,p2eV,p3eV,pfree,pE0],...
    sprintf('Fixed {\\itm}_\\nu^2 = 0 eV^2'),...
    sprintf('Free  {\\itm}_\\nu^2 - \\sigma({\\itm}_\\nu^2) = 0.5  eV^2'),...
    sprintf('Free  {\\itm}_\\nu^2 - \\sigma({\\itm}_\\nu^2) = 1  eV^2'),...
    sprintf('Free  {\\itm}_\\nu^2 - \\sigma({\\itm}_\\nu^2) = 2 eV^2'),...
    sprintf('Free  {\\itm}_\\nu^2 - \\sigma({\\itm}_\\nu^2) = 3 eV^2'),...
    sprintf('Free  {\\itm}_\\nu^2 - unconstrained'));%,...
   % sprintf('Free  {\\itm}_\\nu^2 - unconstrained - \\sigma({\\itE}_0) = 1 eV'));

 
 switch S.RunAnaObj.DataType
     case 'Real'
         leg.Title.String = sprintf('KATRIN KSN1 %.0f%% C.L.',S.ConfLevel);
     case 'Twin'
             leg.Title.String = sprintf('KATRIN MC KSN1 %.0f%% C.L.',S.ConfLevel);
 end
leg.Title.FontWeight = 'normal';
 % sprintf('KATRIN KSN1 %.0f%% C.L. - {\\itm}_\\nu^2 free - \\sigma({\\itm}_\\nu^2) = 1.94 eV^2',S.ConfLevel),... 
legend boxoff;

%% style
PrettyFigureFormat
set(gca,'FontSize',20);
leg.FontSize = 15;
t = get(gca,'title');
t.FontSize = 18;

xlim([6e-03 0.4])
ylim([1 2e3])

title(sprintf('%s - %.0f eV range',leg.Title.String,S.range))
leg.Title.String = '';
%% save
plotname = strrep(S.DefPlotName,'_pull17',sprintf('_ComparemNuPulls_%s.png',S.RunAnaObj.DataType));
print(gcf,plotname,'-dpng','-r300');
fprintf('save plot to %s \n',plotname);
export_fig(gcf,[plotname,'.pdf']);