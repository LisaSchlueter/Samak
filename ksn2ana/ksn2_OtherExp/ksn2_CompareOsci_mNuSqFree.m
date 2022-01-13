% ksn2 calculate chi2 grid search
% compare m2 free, m2 nuisance parameter
%% settings that might change
chi2 = 'chi2CMShape';
DataType = 'Real';%'Twin';%'Real';
nGridSteps = 40;
range = 40;
Mode = 'Compute';
LocalFontSize = 24;
% configure RunAnalysis object
if strcmp(chi2,'chi2Stat')
    NonPoissonScaleFactor = 1;
elseif  strcmp(chi2,'chi2CMShape')
    NonPoissonScaleFactor = 1.112;
end
RunAnaArg = {'RunList','KNM2_Prompt',...
    'DataType',DataType,...
    'fixPar','mNu E0 Norm Bkg',...%free par
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
% configure Sterile analysis object
SterileArg = {'RunAnaObj',A,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
    'nGridSteps',nGridSteps,...
    'SmartGrid','OFF',...
    'RecomputeFlag','OFF',...
    'SysEffect','all',...
    'RandMC','OFF',...
    'range',range,...
    'LoadGridArg',{'mNu4SqTestGrid',5}};

S = SterileAnalysis(SterileArg{:});
S.InterpMode = 'spline';
%%
S.LoadGridFile(S.LoadGridArg{:},'ExtmNu4Sq','ON'); 
S.Interp1Grid('nInter',1e3);
S.FindBestFit;
mNu4Sq_bf = S.mNu4Sq_bf;
sin2T4_bf = S.sin2T4_bf;

PlotPar = S.mNuSq;

%% convert into oscillation parameter space
[~,sin2T4Sq] = S.Convert2Osci;

Ue2Sq = 0.544^2; % from book: particles and nulcei, probably not most up to date
Ue3Sq = 0.151^2; % from book: particles and nulcei, probably not most up to date
Deltam21Sq = 7.37.*1e-05;
Deltam31Sq_NO = 2.56.*1e-03;
Deltam31Sq_IO = -2.54.*1e-03;

Deltam41Sq_NO = S.mNu4Sq-(S.mNuSq-(Ue2Sq.*Deltam21Sq+Ue3Sq.*Deltam31Sq_NO))./(1-S.sin2T4);
Deltam41Sq_IO = S.mNu4Sq-(S.mNuSq-(Ue2Sq.*Deltam21Sq+Ue3Sq.*Deltam31Sq_IO))./(1-S.sin2T4);

%%
close all
GetFigure;
pNO = contour(sin2T4Sq,Deltam41Sq_NO,S.chi2-S.chi2_ref,[5.99 5.99],'LineWidth',2,'Color',rgb('Black'));
hold on;
pIO = contour(sin2T4Sq,Deltam41Sq_IO,S.chi2-S.chi2_ref,[5.99 5.99],':','LineWidth',2,'Color',rgb('Orange'));
pReg = contour(sin2T4Sq,S.mNu4Sq,S.chi2-S.chi2_ref,[5.99 5.99],'-.','LineWidth',2);

PrettyFigureFormat;
set(gca,'XScale','log');
set(gca,'YScale','log');
ylim([0.1 1600]);

xlabel(sprintf('sin^2(2\\theta)'));
ylabel(sprintf('\\Deltam_{41}^2'));
