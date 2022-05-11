% ksn2 calculate chi2 grid search
%% settings that might change
DataType              = 'Real';
range                 = 40;
chi2                  = 'chi2CMShape';
freePar               = 'mNu E0 Norm Bkg';


if contains(freePar,'mNu')
    nGridSteps = 40;
else
    nGridSteps = 30;
end

if strcmp(chi2,'chi2Stat')
    NonPoissonScaleFactor = 1;
elseif  strcmp(chi2,'chi2CMShape')
    NonPoissonScaleFactor = 1.064;
end

RunAnaArg = {'RunList','KNM1',...
    'chi2',chi2,...
    'DataType',DataType,...
    'fixPar',freePar,...
    'NonPoissonScaleFactor',NonPoissonScaleFactor,...
    'SysBudget',200,...
    'minuitOpt','min ; minos',...
    'FSDFlag','KNM2_0p1eV',...
    'ELossFlag','KatrinT2A20',...
    'AngularTFFlag','ON',...
    'SynchrotronFlag','ON',...
    'RadiativeFlag','ON',...
    'DopplerEffectFlag','FSD',...
    'BKG_PtSlope',-2.2*1e-06};


%% configure RunAnalysis object
Real = MultiRunAnalysis(RunAnaArg{:});

% define fit ranges
RangeStandard = Real.GetexclDataStart(40);
qU = round(Real.RunData.qU-18575);
ranges = sort(round(-qU(1:RangeStandard)));


%%
for i=numel(ranges)%1:numel(ranges)
    if i>=2
        % Real = MultiRunAnalysis(RunAnaArg{:});
    end
    Real.exclDataStart = Real.GetexclDataStart(ranges(i));
    % configure Sterile analysis object
    SterileArg = {'RunAnaObj',Real,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
        'nGridSteps',nGridSteps,...
        'SmartGrid','OFF',...
        'RecomputeFlag','OFF',...
        'SysEffect','all',...
        'RandMC','OFF',...
        'range',ranges(i)};
    
    S = SterileAnalysis(SterileArg{:});
    S.GridSearch('ExtmNu4Sq','OFF','mNu4SqTestGrid',2);
    
end



