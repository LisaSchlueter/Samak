% Test SterileAnalysis class
% Lisa, May2020
% plot with neutrino mass as
% nuissance parameter free
% nuissance parameter + 1 eV^2 mnu pull
% nuissance parameter + 2 eV^2 mnu pull
% nuissance parameter + 3 eV^2 mnu pull
%% settings for runanalysis
DataType = 'Real';
%%
RunAnaArg = {'RunList','KNM1',...
    'fixPar','E0 Norm Bkg',...
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
S.RunAnaObj.chi2 = 'chi2Stat';
S.LoadGridFile('CheckSmallerN','ON','CheckLargerN','ON'); % if CheckSmallerN also look for grid with more/less nGridSteps
S.Interp1Grid('RecomputeFlag','ON');% interpolate chi2 map -> nicer appearance of all plots. some
pStat =S.ContourPlot('HoldOn','OFF','Color',rgb('DodgerBlue'),'LineStyle','-');

S.RunAnaObj.chi2 = 'chi2StatNP';
S.LoadGridFile('CheckSmallerN','ON','CheckLargerN','ON'); % if CheckSmallerN also look for grid with more/less nGridSteps
S.Interp1Grid('RecomputeFlag','ON');% interpolate chi2 map -> nicer appearance of all plots. some
pNP = S.ContourPlot('HoldOn','ON','Color',rgb('Orange'),'LineStyle','-.');

% S.RunAnaObj.chi2 = 'chi2CMShape';
% S.LoadGridFile('CheckSmallerN','ON','CheckLargerN','ON'); % if CheckSmallerN also look for grid with more/less nGridSteps
% S.Interp1Grid('RecomputeFlag','ON');% interpolate chi2 map -> nicer appearance of all plots. some
% pCM = S.ContourPlot('HoldOn','ON','Color',rgb('SeaGreen'),'LineStyle',':');

%%
leg = legend([pStat,pNP],...
   ...% sprintf('Fixed {\\itm}_\\nu^2 = 0 eV^2'),...
    sprintf('Stat. only'),...
    sprintf('Stat. and NP background'));
  

 switch S.RunAnaObj.DataType
     case 'Real'
         leg.Title.String = sprintf('KATRIN KSN1 - %.0f eV range at %.0f%% C.L.',S.range,S.ConfLevel);
     case 'Twin'
        leg.Title.String = sprintf('KATRIN MC KSN1 - %.0f eV range at %.0f%% C.L.',S.range,S.ConfLevel);
 end
leg.Title.FontWeight = 'normal';
 % sprintf('KATRIN KSN1 %.0f%% C.L. - {\\itm}_\\nu^2 free - \\sigma({\\itm}_\\nu^2) = 1.94 eV^2',S.ConfLevel),... 
legend boxoff;

%% style
PrettyFigureFormat
set(gca,'FontSize',20);
leg.FontSize = 15;
%t = get(gca,'title');
%t.FontSize = 18;

xlim([6e-03 0.5])
ylim([1 1.7e3])

title('');%sprintf('%.0f eV range',S.range))
%% save
plotname = [S.DefPlotName,'.png'];
print(gcf,plotname,'-dpng','-r300');
fprintf('save plot to %s \n',plotname);

