% calculate chi^2 at some grid point

%% kafit best fits
mNu4Sq = 92.98;
sin2T4 = 0.026;
% result: chi2min = 25.03
%% samak best fits 
mNu4Sq = 97.3;
sin2T4 = 0.027;
% result: chi2min = 24.99
%% something else
mNu4Sq = 95;
sin2T4 = 0.027;
% result: chi2min = 25.005
%% something else
mNu4Sq = 90;
sin2T4 = 0.026;
% result: chi2min = 25.09
%%
freePar  = 'mNu E0 Norm Bkg';
chi2     = 'chi2CMShape';
DataType = 'Real';
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
A.exclDataStart = A.GetexclDataStart(40);
%%
A.ModelObj.SetFitBiasSterile(mNu4Sq,sin2T4);
A.Fit;