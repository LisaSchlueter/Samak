% ksn2 calculate chi2 grid search for all ranges
%% settings that might change
chi2 = 'chi2CMShape';
DataType = 'Real';
range = 40;
freePar = 'E0 Norm Bkg';
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
Real = MultiRunAnalysis(RunAnaArg{:});
%% configure Sterile analysis object


%% define fit ranges
RangeStandard = Real.GetexclDataStart(40);
qU = round(Real.RunData.qU-18574);
ranges = sort(round(-qU(1:RangeStandard)));

if strcmp(DataType,'Real')
    if contains(freePar,'mNu')
        nGridSteps = 40;
        LoadGridArg = {'ExtmNu4Sq','OFF','mNu4SqTestGrid',5} ;
    else
        nGridSteps = 30;
        LoadGridArg = {'ExtmNu4Sq','ON','mNu4SqTestGrid',5} ;
    end
else
    if contains(freePar,'mNu')
        nGridSteps = 40;
    else
        nGridSteps = 30;
    end
    LoadGridArg = {'ExtmNu4Sq','OFF','mNu4SqTestGrid',5} ;
end

   SterileArg = {'RunAnaObj',Real,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
        'SmartGrid','OFF',...
        'RecomputeFlag','OFF',...
        'SysEffect','all',...
        'RandMC','OFF',... 
        'LoadGridArg',LoadGridArg,...
        'nGridSteps',nGridSteps };
    
for i=ranges(end)%numel(ranges):-1:1
    if i<numel(ranges)
        Real = MultiRunAnalysis(RunAnaArg{:}); % re-init for sanity
    end
    Real.exclDataStart = Real.GetexclDataStart(ranges(i));
    % configure Sterile analysis object

    S = SterileAnalysis(SterileArg{:},'range',ranges(i));
    S.GridSearch(S.LoadGridArg{:});
end

