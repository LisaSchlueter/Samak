% ksn1 - study influence of energy-loss parametrizations
% Lisa, May 2020

% work flow:

%% 1. RunAnalysis object
RunAnaArg = {'RunList','KNM1',...
    'fixPar','E0 Norm Bkg',...
    'DataType','Real',...
    'FSDFlag','SibilleFull',...
    'ELossFlag','KatrinT2',...
    'AnaFlag','StackPixel',...
    'chi2','chi2Stat',...
    'ROIFlag','Default',...
    'SynchrotronFlag','ON',...
    'AngularTFFlag','OFF',...
    'ISCSFlag','Edep',...
    'TwinBias_Q',18573.73,...
    'SysBudget',29,...
    'pullFlag',99,...
    'NonPoissonScaleFactor',1};

R = MultiRunAnalysis(RunAnaArg{:});
R.chi2 = 'chi2CMShape';
%% 2. SterileAnalysis class
SterileArg = {'RunAnaObj',R,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
    'nGridSteps',50,...
    'SmartGrid','OFF',...
    'RecomputeFlag','OFF',...
    'SysEffect','all',...
    'RandMC','OFF',...
    'range',95};

S = SterileAnalysis(SterileArg{:});
S.RunAnaObj.SysBudget= 29;
S.RunAnaObj.ELossFlag = 'KatrinT2';
S.RunAnaObj.AngularTFFlag = 'OFF';
%% plot
Mode = 'RhoDerr';
S.range = 95;
CL = 95;

% new: e-loss, different uncertainty etc.
switch Mode
    case 'Eloss'
        S.RunAnaObj.ELossFlag = 'KatrinT2A20';
        legStrNew = 'KNM2 KATRIN energy-loss function (April 2020)';
        legStrOld = 'KNM1 KATRIN energy-loss function';
    case 'RhoDerr'
        S.RunAnaObj.SysBudget= 254;
        legStrNew = sprintf('\\Delta\\rhod\\sigma = 2.00%%');
        legStrOld = sprintf('\\Delta\\rhod\\sigma = 0.85%%');
    case 'Plasmaerr'
        S.RunAnaObj.SysBudget = 299;
        legStrNew = sprintf('\\Delta\\sigma^2 = 0.0021 eV^2, \\Delta\\epsilon_{loss} = 0.04 eV');
        legStrOld = sprintf('Broadening + e-loss shift neglected');
    case 'AngTF'
        S.RunAnaObj.AngularTFFlag = 'ON';
        legStrNew = sprintf('With non-isotropic transmission function');
        legStrOld = sprintf('Without non-isotropic transmission function');
end
S.LoadGridFile('CheckSmallerN','ON','CheckLargerN','ON'); % if CheckSmallerN also look for grid with more/less nGridSteps
S.Interp1Grid('RecomputeFlag','ON');% interpolate chi2 map -> nicer appearance of all plots. some
pnew = S.ContourPlot('BestFit','ON','SavePlot','OFF','CL',CL,'HoldOn','OFF');

% default case:
S.RunAnaObj.SysBudget= 29;
S.RunAnaObj.ELossFlag = 'KatrinT2';
S.RunAnaObj.AngularTFFlag = 'OFF';
S.LoadGridFile('CheckSmallerN','ON','CheckLargerN','ON'); % if CheckSmallerN also look for grid with more/less nGridSteps
S.InterpMode = 'lin'; % waring: if contour is closed, spline interp sometimes sensitive to artefacts! Switch to "lin" in this case
S.Interp1Grid('RecomputeFlag','ON');% interpolate chi2 map -> nicer appearance of all plots. some
pold = S.ContourPlot('BestFit','ON','SavePlot','OFF','CL',CL,'Color',rgb('Orange'),'LineStyle','-.','HoldOn','ON');

xlim([1e-03 0.5])
ylim([1 1e4])
leg = legend([pnew,pold],legStrNew,legStrOld,...
    'EdgeColor',rgb('Silver'),'Location','southwest');

title(sprintf('%s - %.0f%% C.L.',S.GetPlotTitle,CL),'FontWeight','normal','FontSize',get(gca,'FontSize'));

plotname = sprintf('%s_%sComparison_%.2gCL.png',S.DefPlotName,Mode,CL);
print(gcf,plotname,'-dpng','-r450');
fprintf('save plot to %s \n',plotname);