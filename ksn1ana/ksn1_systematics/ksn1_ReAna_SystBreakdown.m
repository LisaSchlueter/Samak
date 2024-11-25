% ksn2 calculate chi2 grid search
%% settings that might change
nGridSteps            = 30;
DataType              = 'Twin';
range                 = 40;
chi2                  = 'chi2Stat';
freePar               = 'E0 Norm Bkg';
%% configure RunAnalysis object

A = MultiRunAnalysis('RunList','KNM1',...
    'chi2',chi2,...
    'DataType',DataType,...
    'fixPar',freePar,...
    'NonPoissonScaleFactor',1,...
    'SysBudget',200,...
    'minuitOpt','min ; minos',...
    'FSDFlag','KNM2_0p1eV',...
    'ELossFlag','KatrinT2A20',...
    'AngularTFFlag','ON',...
    'SynchrotronFlag','ON',...
    'RadiativeFlag','ON',...
    'DopplerEffectFlag','FSD',...
    'BKG_PtSlope',-2.2*1e-06);
A.exclDataStart = A.GetexclDataStart(range);
%% configure Sterile analysis object
SterileArg = {'RunAnaObj',A,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
    'nGridSteps',nGridSteps,...
    'SmartGrid','OFF',...
    'RecomputeFlag','OFF',...
    'SysEffect','all',...
    'RandMC','OFF',...
    'range',range};
S = SterileAnalysis(SterileArg{:});
S.LoadGridArg = {'ExtmNu4Sq','OFF','mNu4SqTestGrid',2};

S.GridSearch('ExtmNu4Sq','OFF','mNu4SqTestGrid',2);

return
SysEffectsAll   = {'Stat','BkgPT','NP','FSD', 'RF_BF','Bkg','TASR','RF_EL',...
    'RF_RX','FPDeff','Stack','TCoff_OTHER','all'};

for i=1:numel(SysEffectsAll)
    S.SysEffect = SysEffectsAll{i};
    
    if strcmp(SysEffectsAll{i},'NP') || strcmp(SysEffectsAll{i},'Stat')
        S.RunAnaObj.chi2 = 'chi2Stat';
        if strcmp(SysEffectsAll{i},'NP')
            S.RunAnaObj.NonPoissonScaleFactor = 1.064;
        end
        S.GridSearch(S.LoadGridArg{:});
        S.RunAnaObj.chi2 = 'chi2CMShape';
        S.RunAnaObj.NonPoissonScaleFactor = 1;
    elseif strcmp(SysEffectsAll{i},'all')
        S.RunAnaObj.NonPoissonScaleFactor = 1.064;
        S.GridSearch(S.LoadGridArg{:});
        S.RunAnaObj.NonPoissonScaleFactor = 1;
    else
        S.GridSearch(S.LoadGridArg{:});
    end
end

