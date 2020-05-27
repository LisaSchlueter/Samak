% ksn1 - study influence of some systematics/modelling
% options: 
% new (KNM2 T2) energy-loss vs. old (KNM1 T2) parametrizations
% plasma uncertainty vs. no plasma uncertainty
% enhance column density uncertainty: 0.85% -->1.5% or 2%
% enhance FSD onset uncertainty: 1% --> 1.5%
% enhance FSD excited states uncertainty: 18% --> 50%
% non-isotropic transmission function 

% Lisa, May 2020

%% 1. RunAnalysis object
RunAnaArg = {'RunList','KNM1',...
    'fixPar','mNu E0 Norm Bkg',...
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
    'NonPoissonScaleFactor',1.064};

A = MultiRunAnalysis(RunAnaArg{:});
A.chi2 = 'chi2Stat';
%% 2. SterileAnalysis class
SterileArg = {'RunAnaObj',R,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
    'nGridSteps',50,...
    'SmartGrid','OFF',...
    'RecomputeFlag','OFF',...
    'SysEffect','all',...
    'RandMC','OFF',...
    'range',95};

S = SterileAnalysis(SterileArg{:});
S.RunAnaObj.SysBudget= 24;
S.RunAnaObj.ELossFlag = 'KatrinT2';
S.RunAnaObj.AngularTFFlag = 'OFF';
%% plot
CL = 96.46;
S.InterpMode = 'lin'; % waring: if contour is closed, spline interp sometimes sensitive to artefacts! Switch to "lin" in this case
S.range = 95;

Mode = 'AngTFEloss';

switch Mode
    case 'Eloss' % new: e-loss, different uncertainty etc.
        S.RunAnaObj.ELossFlag = 'KatrinT2A20';
        legStrNew = 'KNM2 KATRIN energy-loss function (April 2020)';
        legStrOld = 'KNM1 KATRIN energy-loss function';
    case 'RhoDerr2'
        S.RunAnaObj.SysBudget= 254;
        legStrNew = sprintf('\\Delta\\rhod\\sigma = 2.00%%');
        legStrOld = sprintf('\\Delta\\rhod\\sigma = 0.85%%');
    case 'RhoDerr1p5'
        S.RunAnaObj.SysBudget= 256;
        legStrNew = sprintf('\\Delta\\rhod\\sigma = 1.50%%');
        legStrOld = sprintf('\\Delta\\rhod\\sigma = 0.85%%');
    case 'Plasmaerr'
        S.RunAnaObj.SysBudget = 299;
        legStrNew = sprintf('\\Delta\\sigma^2 = 0.0021 eV^2, \\Delta\\epsilon_{loss} = 0.04 eV');
        legStrOld = sprintf('Broadening + e-loss shift neglected');
    case 'AngTF'
        S.RunAnaObj.AngularTFFlag = 'ON';
        legStrNew = sprintf('With non-isotropic transmission function');
        legStrOld = sprintf('Without non-isotropic transmission function');
    case 'AngTFEloss'
        S.RunAnaObj.AngularTFFlag = 'ON';
        S.RunAnaObj.ELossFlag = 'KatrinT2A20';
        legStrNew = sprintf('With non-isotropic transmission function + KNM2 energy-loss');
        legStrOld = sprintf('Without non-isotropic transmission function + KNM1 energy-loss');
    case 'FSDonset'
        S.RunAnaObj.SysBudget = 256;
        legStrNew = sprintf('\\Delta FSD onset = 1.5%%');
        legStrOld = sprintf('\\Delta FSD onset = 1.0%%');
    case 'FSDexcitedStates50'
        S.RunAnaObj.SysBudget = 257;
        legStrNew = sprintf('\\Delta FSD excited states (bin-to-bin) = 50%%');
        legStrOld = sprintf('\\Delta FSD excited states (bin-to-bin) = 18%%');
    case 'FSDexcitedStates45'
        S.RunAnaObj.SysBudget = 258;
        legStrNew = sprintf('\\Delta FSD excited states (bin-to-bin) = 45%%');
        legStrOld = sprintf('\\Delta FSD excited states (bin-to-bin) = 18%%');
    case 'FSDexcitedStates40'
        S.RunAnaObj.SysBudget = 253;
        legStrNew = sprintf('\\Delta FSD excited states (bin-to-bin) = 40%%');
        legStrOld = sprintf('\\Delta FSD excited states (bin-to-bin) = 18%%');
end
S.LoadGridFile('CheckSmallerN','ON','CheckLargerN','ON'); % if CheckSmallerN also look for grid with more/less nGridSteps
S.Interp1Grid('RecomputeFlag','ON');% interpolate chi2 map -> nicer appearance of all plots. some
pnew = S.ContourPlot('BestFit','ON','SavePlot','OFF','CL',CL,'HoldOn','OFF');

% default ksn1 contour:
S.RunAnaObj.SysBudget= 24;
S.RunAnaObj.ELossFlag = 'KatrinT2';
S.RunAnaObj.AngularTFFlag = 'OFF';
S.LoadGridFile('CheckSmallerN','ON','CheckLargerN','ON'); % if CheckSmallerN also look for grid with more/less nGridSteps
S.Interp1Grid('RecomputeFlag','ON');% interpolate chi2 map -> nicer appearance of all plots. some
pold = S.ContourPlot('BestFit','ON','SavePlot','OFF','CL',CL,'Color',rgb('Orange'),'LineStyle','-.','HoldOn','ON');

xlim([1e-03 0.5])
ylim([1 1e4])
leg = legend([pnew,pold],legStrNew,legStrOld,...
    'EdgeColor',rgb('Silver'),'Location','southwest');

%title(sprintf('%s - %.1f%% C.L.',S.GetPlotTitle,CL),'FontWeight','normal','FontSize',get(gca,'FontSize'));
title(sprintf('%s - \\Delta\\chi^2 = %.2f',S.GetPlotTitle,GetDeltaChi2(CL,2)),'FontWeight','normal','FontSize',get(gca,'FontSize'));

plotname = sprintf('%s_%sComparison_%.2gCL.png',S.DefPlotName,Mode,CL);
print(gcf,plotname,'-dpng','-r450');
fprintf('save plot to %s \n',plotname);
export_fig(gcf,strrep(plotname,'.png','.pdf'));