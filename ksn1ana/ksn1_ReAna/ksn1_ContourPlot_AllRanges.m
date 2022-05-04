% ksn2 calculate chi2 grid search
%% settings that might change
nGridSteps            = 30;
DataType              = 'Real';
range                 = 40;
chi2                  = 'chi2CMShape';
freePar               = 'E0 Norm Bkg';



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

%% define fit ranges
RangeStandard = Real.GetexclDataStart(40);
qU = round(Real.RunData.qU-18575);
ranges = sort(round(-qU(1:RangeStandard)));

%ranges = [ranges_All(1:9);ranges_All(12:13)];
  SterileArg = {'RunAnaObj',Real,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
        'nGridSteps',nGridSteps,...
        'SmartGrid','OFF',...
        'RecomputeFlag','OFF',...
        'SysEffect','all',...
        'RandMC','OFF',...
        'range',40,...
        'LoadGridArg',{'ExtmNu4Sq','OFF','mNu4SqTestGrid',2},...
        'ConfLevel',0.99};
     S = SterileAnalysis(SterileArg{:});
    
 %
for i=1:numel(ranges)
    Real.exclDataStart = Real.GetexclDataStart(ranges(i));
    % configure Sterile analysis object
  
    S.range = ranges(i);
    S.LoadGridFile(S.LoadGridArg{:});
   
   S.Interp1Grid;
   if i==1
       S.ContourPlot('BestFit','ON');
   else
       S.ContourPlot('BestFit','ON','HoldOn','ON','Color',S.PlotColors{i},'LineStyle',S.PlotLines{i});
   end
    
    
end

%%

Real.exclDataStart = Real.GetexclDataStart(ranges(end));
% configure Sterile analysis object
S.ConfLevel = 99;
S.range = ranges(end);
%%
S.RunAnaObj.FSDFlag = 'KNM2_0p1eV';
S.LoadGridFile(S.LoadGridArg{:});
S.InterpMode = 'spline';
S.Interp1Grid;
S.ContourPlot('BestFit','ON');

S.RunAnaObj.FSDFlag = 'KNM2';
S.LoadGridFile(S.LoadGridArg{:});
S.InterpMode = 'spline';
S.Interp1Grid;
S.ContourPlot('BestFit','ON','HoldOn','ON','Color',rgb('Orange'),'LineStyle','-.','MarkerStyle','o');

 
 