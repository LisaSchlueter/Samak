% save my contour to .dat file for other fitters to read

% ksn2 calculate chi2 grid search
%% settings that might change
chi2 = 'chi2Stat';%CMShape';
DataType = 'Twin';
nGridSteps = 40;
range = 40;
freePar = 'mNu E0 Norm Bkg';
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
SterileArg = {'RunAnaObj',A,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
    'nGridSteps',nGridSteps,...
    'SmartGrid','OFF',...
    'RecomputeFlag','OFF',...
    'SysEffect','all',...
    'RandMC','OFF',...
    'range',range};
S = SterileAnalysis(SterileArg{:});
%%
close all
S.LoadGridFile('mNu4SqTestGrid',5);
S.InterpMode = 'spline';
S.Interp1Grid('Maxm4Sq',38.2^2);
S.ContourPlot; %close;

% S.InterpMode = 'lin';
% S.LoadGridFile;
% S.Interp1Grid('Maxm4Sq',38^2);
% S.ContourPlot('HoldOn','ON','Color',rgb('Orange'),'LineStyle','--'); %close;
% 
% mNu4Sq   = linspace(min(S.mNu4Sq_contour),max(S.mNu4Sq_contour),1e4);
% sin2T4   = interp1(S.mNu4Sq_contour,S.sin2T4_contour,mNu4Sq,'spline');
% 
% hold on;
% plot(sin2T4,mNu4Sq)
% S.sin2T4_contour = sin2T4;
% S.mNu4Sq_contour = mNu4Sq;
%% export contour

savedir = sprintf('%sksn2ana/ksn2_FitterCrossCheck/results/',getenv('SamakPath'));
MakeDir(savedir);
savename =  sprintf('%sKSN2_contour_Samak%s_%s_40eV_%s',...
    savedir,DataType,strrep(freePar,' ',''),strrep(chi2,'chi2',''));
Write2Txt('filename',savename,...
    'Format','dat','variable',[S.sin2T4_contour;S.mNu4Sq_contour],'nCol',2,'variableName','sinT4Sq m4Sq');
%% test
d = importdata([savename,'.dat']);
plot(d.data(:,1),d.data(:,2));
set(gca,'XScale','log');
set(gca,'YScale','log');
xlim([5e-03,0.5]);
ylim([1 40^2]);

