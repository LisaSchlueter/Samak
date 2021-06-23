% ksn2 chi2 grid - non physical parameter space 4 quadrants
% ksn2 plot chi2 grid search
% 4 quadrants
%% settings that might change
chi2 = 'chi2CMShape';
DataType = 'Real';
nGridSteps = 40;
range = 40;

Chi2Crit_H0 = 6.31;%6.19; % +- 0.12 (from Martin Slezak, H0, 5000 contours)
Chi2CritErr_H0 = 0.30;
Chi2Crit_H1 = 6.69;
Chi2CritErr_H1 = 0.27;

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
    'FSDFlag','KNM2_0p5eV',...
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
SterileArg = {... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
    'nGridSteps',nGridSteps,...
    'SmartGrid','OFF',...
    'RecomputeFlag','OFF',...
    'SysEffect','all',...
    'RandMC','OFF',...
    'range',range,...
    'LoadGridArg',{'ExtmNu4Sq','ON','mNu4SqTestGrid',5},...
    'InterpMode','spline'};
S = SterileAnalysis('RunAnaObj',A,SterileArg{:});
%%
S.LoadGridFile(S.LoadGridArg{:});
S.Interp1Grid;

% Null hypothesis
Chi2Crit_H0 = 6.31;%6.19; % +- 0.12 (from Martin Slezak, H0, 5000 contours)
Chi2CritErr_H0 = 0.30;
myCL_H0 = GetCL(Chi2Crit_H0);
p2 = S.ContourPlot('CL',myCL_H0,'HoldOn','ON','Color',rgb('Skyblue'),'LineStyle','-.');
mNu4Sq_contour_H0 = S.mNu4Sq_contour;
sin2T4_contour_H0 = S.sin2T4_contour;

myCL_H01 = GetCL(Chi2Crit_H0+Chi2CritErr_H0);
p41 = S.ContourPlot('CL',myCL_H01,'HoldOn','ON','Color',rgb('LightSalmon'),'LineStyle',':'); 
mNu4Sq_contour_H0_Up = S.mNu4Sq_contour;
sin2T4_contour_H0_Up = interp1(mNu4Sq_contour_H0_Up,S.sin2T4_contour,mNu4Sq_contour_H0,'spline');

myCL_H02 = GetCL(Chi2Crit_H0-Chi2CritErr_H0);
p51 = S.ContourPlot('CL',myCL_H02,'HoldOn','ON','Color',rgb('LightSalmon'),'LineStyle',':'); 
mNu4Sq_contour_H0_Down = S.mNu4Sq_contour;
sin2T4_contour_H0_Down = interp1(mNu4Sq_contour_H0_Down,S.sin2T4_contour,mNu4Sq_contour_H0,'spline');

sin2T4_contourErr_H0 = [sin2T4_contour_H0_Up-sin2T4_contour_H0;sin2T4_contour_H0-sin2T4_contour_H0_Down];
close all;
%% Sterile hypothesis
Chi2Crit_H1 = 6.69;
Chi2CritErr_H1 = 0.27;
myCL_H1 = GetCL(Chi2Crit_H1);
p3 = S.ContourPlot('CL',myCL_H1,'HoldOn','ON','Color',rgb('LightSalmon'),'LineStyle',':');
mNu4Sq_contour_H1 = S.mNu4Sq_contour;
sin2T4_contour_H1 = S.sin2T4_contour;

myCL_H11 = GetCL(Chi2Crit_H1+Chi2CritErr_H1);
p4 = S.ContourPlot('CL',myCL_H11,'HoldOn','ON','Color',rgb('LightSalmon'),'LineStyle',':'); 
mNu4Sq_contour_H1_Up = S.mNu4Sq_contour;
sin2T4_contour_H1_Up = interp1(mNu4Sq_contour_H1_Up,S.sin2T4_contour,mNu4Sq_contour_H1,'spline');

myCL_H12 = GetCL(Chi2Crit_H1-Chi2CritErr_H1);
p5 = S.ContourPlot('CL',myCL_H12,'HoldOn','ON','Color',rgb('LightSalmon'),'LineStyle',':'); 
mNu4Sq_contour_H1_Down = S.mNu4Sq_contour;
sin2T4_contour_H1_Down = interp1(mNu4Sq_contour_H1_Down,S.sin2T4_contour,mNu4Sq_contour_H1,'spline');

sin2T4_contourErr_H1 = [sin2T4_contour_H1_Up-sin2T4_contour_H1;sin2T4_contour_H1-sin2T4_contour_H1_Down];
close all;

%% plot
PlotMode = 'Bounded'; %Boundedline
GetFigure;
switch PlotMode
    case 'Plt'  
        p1 = S.ContourPlot('HoldOn','OFF','Color',rgb('Black'),'CL',0.95);
        hold on;
        pH0 = plot(sin2T4_contour_H0,mNu4Sq_contour_H0,'Color',rgb('Orange'),'LineWidth',2,'LineStyle','-.');
        pH1 =  plot(sin2T4_contour_H1,mNu4Sq_contour_H1,'Color',rgb('FireBrick'),'LineWidth',2,'LineStyle',':');
         legHandles = [p1,pH0,pH1];
    case 'Bounded'
        [pH1,a] = boundedline(sin2T4_contour_H1,mNu4Sq_contour_H1,sin2T4_contourErr_H1','orientation','horiz');
        pH1.Color = rgb('FireBrick'); a.FaceColor = rgb('FireBrick');  a.FaceAlpha  =0.4;
        pH1.LineWidth = 1.5;pH1.LineStyle = ':';
        hold on;
        [pH0,a2] = boundedline(sin2T4_contour_H0,mNu4Sq_contour_H0,sin2T4_contourErr_H0','orientation','horiz');
        pH0.Color = rgb('Orange'); a2.FaceColor = rgb('Orange'); a2.FaceAlpha  =0.4;
        pH0.LineWidth = 1.5;pH0.LineStyle = '-.';
        hold on;
        p1 = S.ContourPlot('HoldOn','ON','Color',rgb('Black'),'CL',0.95);
        p1.LineWidth = 1.5;
        legHandles = [p1,a2,a];
end

xlim([5e-03, 0.5]);
ylim([0.5, 1600]);
title('');
leg = legend(legHandles,sprintf('\\chi^2_{crit.} = 5.99 (Wilk`s theorem)'),...
    sprintf('\\chi^2_{crit.} = %.2f \\pm %.2f (MC simulation: {\\itm}_4^2 = 0 eV^2 , |{\\itU}_{e4}|^2 = 0)',Chi2Crit_H0,Chi2CritErr_H0),...
      sprintf('\\chi^2_{crit.} = %.2f \\pm %.2f (MC simulation: {\\itm}_4^2 = 92.7 eV^2 , |{\\itU}_{e4}|^2 = 0.024)',Chi2Crit_H1,Chi2CritErr_H1));
PrettyLegendFormat(leg);
legend box off

savedir = [getenv('SamakPath'),'ksn2ana/ksn2_WilksTheorem/plots/'];
MakeDir(savedir);
pltname = sprintf('%sksn2_WT_ContourShift_Chi2Crit%.2f.pdf',savedir,Chi2Crit_H0);
export_fig(pltname);
fprintf('save plot to %s \n',pltname)
%%
function cl = GetCL(chi2crit)
x = linspace(80,99.99,1e3);
chi2 = GetDeltaChi2(x,2);
cl = interp1(chi2,x,chi2crit);
end
