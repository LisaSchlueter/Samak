% ksn2 calculate chi2 grid search
% add file with small mNu4Sq
%% settings that might change
chi2 = 'chi2Stat';
DataType = 'Real';
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
%% calculate both: normal grid + additional small m4^2 grid
S.GridSearch('AddSmallNu4Sq','OFF');
S.GridSearch('AddSmallNu4Sq','ON');

%% load both
S.LoadGridFile('AddSmallNu4Sq','ON');
chi2_s   = S.chi2(2:end,:); % 10x25
mNu4Sq_s = S.mNu4Sq(2:end,:);
sin2T4_s = S.sin2T4(2:end,:);
mNuSq_s  = S.mNuSq(2:end,:);
E0_s     = S.E0(2:end,:);

S.GridPlot('Contour','OFF','BestFit','OFF','SavePlot','OFF','CL',95)

%%
S.LoadGridFile('AddSmallNu4Sq','OFF');
chi2_l   = S.chi2;  % 25 x 25
mNu4Sq_l = S.mNu4Sq;
sin2T4_l = S.sin2T4;
mNuSq_l  = S.mNuSq;
E0_l     = S.E0;

S.Interp1Grid('RecomputeFlag','ON');

S.GridPlot('Contour','OFF','BestFit','OFF','SavePlot','OFF','CL',95)

%%
S.mNu4Sq = [mNu4Sq_s;mNu4Sq_l];
S.sin2T4 = [sin2T4_s;sin2T4_l];
S.chi2 = [chi2_s;chi2_l]';
S.mNuSq = [mNuSq_s;mNuSq_l]';
S.E0 = [E0_s;E0_l]';
S.InterpMode = 'lin';
S.Interp1Grid('RecomputeFlag','ON');
S.chi2 = S.chi2';
S.GridPlot('Contour','OFF','BestFit','OFF','SavePlot','OFF','CL',95)

% if strcmp(A.DataType,'Real')
%     S.InterpMode = 'lin';
%     BF = 'ON';
% else
%     S.InterpMode = 'spline';
%    BF = 'OFF';
% end
% S.Interp1Grid('RecomputeFlag','ON');
%
%S.ContourPlot('BestFit','OFF','CL',95)
% S.PlotStatandSys('SavePlot','png')
%S.PlotmNuSqOverview('PullmNuSq','OFF','SavePlot','png')