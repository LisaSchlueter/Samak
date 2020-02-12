F=FakeRunsGenerator('NumberOfRuns',90,'FakeRunPath','../../tritium-fakedata/mat');
F.DrawFakeRunsTBDIS_Uniform;
F.SaveFakeRuns;
stairs(F.WGTS_MolFrac_TT_LuT); PrettyFigureFormat
F.SimFakeObj_LuT{5}.PlotTBDIS('CPS','ON')

return;

% w=RunAnalysis('DataType','Fake','FakeRunType','Fake1','RunNr',240);
% w.Fit; 
% w.PlotFit;
% w.PlotDataModel_KNM1

% M=MultiRunAnalysis('DataType','Fake','FakeRunType','Fake1','RunList','KNM1_30d','exclDataStart',1,'fixPar','5 6','chi2','chi2Stat','Debug','ON','ELossFlag','Abdurashitov');
% M.FitRunList; 

M=MultiRunAnalysis('DataType','Fake','FakeRunType','Fake5y','RunList','KNM5y','exclDataStart',1,'fixPar','5 6','chi2','chi2Stat','Debug','ON','ELossFlag','Abdurashitov');
M.FitRunList; 

% M.StackRuns('CutOnSC','OFF','SCsigma',3,'CutOnFitSingleRuns','OFF');
% M.SimulateStackRuns();
% M.StackRuns('CutOnSC','ON','SCsigma',3,'CutOnFitSingleRuns','ON','Fitsigma',8);

defaultEffects = struct(...
                'RF_EL','OFF',...  % Response Function(RF) EnergyLoss
                'RF_BF','OFF',...  % RF B-Fields
                'RF_RX','ON',...  % Column Density, inel cross ection
                'FSD','ON',...
                'TASR','ON',...
                'TCoff_RAD','ON',...
                'TCoff_OTHER','OFF',...
                'DOPoff','OFF');
M.ComputeCM('RecomputeFlag','OFF','Stack','OFF','WGTS_TASR_RelErr',1e-3,'SysEffects',defaultEffects,'BkgCM','ON');

M.Fit; M.PlotFit('Mode','Rate');
