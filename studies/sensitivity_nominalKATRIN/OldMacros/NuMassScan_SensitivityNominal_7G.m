TimeSec = 3*365*24*60*60;
MACE_Ba_T = 7*1e-04;
WGTS_B_T = 0.7*3.6;
Q_i = 18575;
range = [30,45,60];
Scan = 'ON';
PlotFit = 'OFF'; %only used when Scan ON
SysBudget = '03';

FPD_MeanEff = 0.9;
Arg = {'TimeSec',TimeSec,'Q_i',Q_i,'range',range(1),...
            'MACE_Ba_T',MACE_Ba_T,'WGTS_B_T',WGTS_B_T,...
            'Scan',Scan,'SysBudget',SysBudget,...
            'TD','MTDcreator',...
            'FPD_MeanEff',FPD_MeanEff};        
S = SensitivityStudy(Arg{:});
%
%for s = 1:numel(SysBudget)
for r=1:numel(range)
LoopArg = {'TimeSec',TimeSec,'Q_i',Q_i,'range',range(r),...
    'MACE_Ba_T',MACE_Ba_T,'WGTS_B_T',WGTS_B_T,...
    'Scan',Scan,'PlotFit',PlotFit,'SysBudget',SysBudget,...
    'RecomputeFlag',RecomputeFlag,'saveResults','ON','Anchor6G',Anchor6G};

% Do stat + sys (one by one)
     NuMassScan_SensitivityNominal_Systematics_Loop(LoopArg{:});
     
% Do stat + N+1
[mNu_Fit, parScan, errScan, chi2minScan, dof ,mNu90,mNumin] = ...
    NuMassScan_SensitivityNominal_Nplus1Systematics_Loop(LoopArg{:});
end
%end
%% Bar Plot
SysBudget  = '03';
NuMassScan_Plot_SensitivityNominal_Nplus1Systematics('SingleSys','ON','SysBudget',SysBudget,...
   'MACE_Ba_T',MACE_Ba_T,'WGTS_B_T',WGTS_B_T,'range',range,'TimeSec',TimeSec,...
    'SingleSys','OFF','TDRSys','OFF');
%% RF Breakdown
SysBudget  = '03';
NuMassScan_PlotMultiBar_RFBreakdown('SingleSys','ON','SysBudget',SysBudget,...
   'MACE_Ba_T',MACE_Ba_T,'WGTS_B_T',WGTS_B_T,'range',range,'TimeSec',TimeSec,...
    'SingleSys','ON');


