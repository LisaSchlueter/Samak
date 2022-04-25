% ksn2 extract systematic contribution
% raster scan
% compare variances
%% settings that might change
chi2       = 'chi2CMShape';
DataType   = 'Twin';
nGridSteps = 30;
range      = 40;
CL         =  chi2cdf(1,1);

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
if strcmp(DataType,'Real')
    LoadGridArg = {'mNu4SqTestGrid',5,'ExtmNu4Sq','ON'};
else
    LoadGridArg = {'mNu4SqTestGrid',5};
end

SterileArg = {'RunAnaObj',A,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
    'nGridSteps',nGridSteps,...
    'SmartGrid','OFF',...
    'RecomputeFlag','OFF',...
    'SysEffect','all',...
    'RandMC','OFF',...
    'range',range,...
    'LoadGridArg',LoadGridArg,...
    'ConfLevel',CL};

%%
S = SterileAnalysis(SterileArg{:});
%%
[Ratio,StatDomFraction,mNu4Sq,sin2t4_Stat,sin2t4_Tot,sin2t4_Sys] = S.StatOverSysKsn2('RasterScan','ON');

savedir = [getenv('SamakPath'),'ksn2ana/ksn2_Systematics/results/'];
savename = sprintf('%sksn2_StatOverSyst_%s.mat',savedir,DataType);

if CL~=0.95
    savename = strrep(savename,'.mat',sprintf('%.2gCL.mat',CL));
end

if ~exist(savename,'file')
    save(savename,'mNu4Sq','sin2t4_Stat','sin2t4_Tot','sin2t4_Sys');
end
%%
 close all;
 GetFigure;
hold on;
plot(mNu4Sq,Ratio,'LineWidth',3);
ylabel('Variance ratio');
ylabel(sprintf('\\sigma_{syst.}^2 / \\sigma_{total}^2'))
xlabel(sprintf('{\\itm}_4^2 (eV^2)'));
set(gca,'XScale','log');
PrettyFigureFormat('FontSize',22);


