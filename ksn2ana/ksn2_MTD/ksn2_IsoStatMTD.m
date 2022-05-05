% look at sensitivity contour from KNM-2 like simulation with flat MTD

FakeInitFile = @ref_KNM2_KATRIN_IsoStatMTD;
FakeInitFile = @ref_KNM2_KATRIN_IsoStatMTD_Bkg0mcps;
range = 40;
freePar = 'E0 Norm Bkg';
nGridSteps = 25;

% tritium run model
F = RunAnalysis('RunNr',1,...
    'DataType','Fake',...
    'FakeInitFile',FakeInitFile,...
    'fixPar',freePar,...
    'SysBudget',40,...
    'AnaFlag','StackPixel',...
    'chi2','chi2Stat',...
    'FSDFlag','Sibille0p5eV',...
    'ELossFlag','KatrinT2A20',...
    'DopplerEffectFlag','FSD',...
    'RadiativeFlag','ON',...
    'minuitOpt','min ; minos');

F.exclDataStart = F.GetexclDataStart(range);


%% configure Sterile analysis object
SterileArg = {'RunAnaObj',F,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
    'nGridSteps',nGridSteps,...
    'SmartGrid','OFF',...
    'RecomputeFlag','OFF',...
    'SysEffect','all',...
    'RandMC','OFF',...
    'range',range};
S = SterileAnalysis(SterileArg{:});
S.GridSearch('mNu4SqTestGrid',2);
return
%%
S.LoadGridFile
S.Interp1Grid;
S.ContourPlot