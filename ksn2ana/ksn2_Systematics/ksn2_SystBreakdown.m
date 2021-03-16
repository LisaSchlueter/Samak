% ksn2 calculate chi2 grid search
% compare m2 free, m2 nuisance parameter
%% settings that might change
chi2 = 'chi2Stat';
DataType = 'Twin';
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

%%
S = SterileAnalysis(SterileArg{:});
S.RunAnaObj.chi2 = 'chi2CMShape';
S.RunAnaObj.NonPoissonScaleFactor = 1;
SysEffectsAll   = {'Stat','FSD','BkgPT','Bkg','NP'};%,'NP','LongPlasma','TASR','RF_EL','RF_BF','RF_RX','FPDeff','Stack','TCoff_OTHER'}; %Bkg has to be last

HoldOn = 'OFF';
pHandle = cell(numel(SysEffectsAll),1);
LineStyle = {'-','-.',':','--','-','-.',':','--','-','-.',':','--'};
Colors = colormap('winter');

for i=1:numel(SysEffectsAll)
    S.SysEffect = SysEffectsAll{i};
    
    if strcmp(SysEffectsAll{i},'NP') || strcmp(SysEffectsAll{i},'Stat')
        S.RunAnaObj.chi2 = 'chi2Stat';
        if strcmp(SysEffectsAll{i},'NP')
            S.RunAnaObj.NonPoissonScaleFactor = 1.112;
        end
        S.GridSearch;
        
        S.LoadGridFile;
        S.Interp1Grid;
        
        S.RunAnaObj.chi2 = 'chi2CMShape';
        S.RunAnaObj.NonPoissonScaleFactor = 1;

    else
        S.GridSearch;
        
        S.LoadGridFile;
        S.Interp1Grid;
    end

   pHandle{i} =  S.ContourPlot('HoldOn',HoldOn,'SavePlot','OFF','Color',Colors(i*floor(256/numel(SysEffectsAll)),:),'LineStyle',LineStyle{i});
    HoldOn = 'ON';
end

%%
xlim([7e-03,0.5])
SysEffectsAll   = {'FSD','BkgPT','Bkg','NP'};%,'NP','LongPlasma','TASR','RF_EL','RF_BF','RF_RX','FPDeff','Stack','TCoff_OTHER'}; %Bkg has to be last

SysEffectLeg    = {'Stat. only';'Final-state distribution';'Background PT';'Background slope'; 'Non-Poisson background'};
leg = legend([pHandle{:}],SysEffectLeg);
PrettyLegendFormat(leg);
leg.Title.String = 'Stat. + 1 syst.';
leg.Title.FontWeight ='normal';

%                         'Tritium activity fluctuations';'Energy-loss function';...
%                         'Magnetic fields';'Source scattering';'HV fluctuations';...
%                         'Long. source potential';...
%                         'Theoretical corrections';...
%                         'Detector efficiency';};
%%

            

