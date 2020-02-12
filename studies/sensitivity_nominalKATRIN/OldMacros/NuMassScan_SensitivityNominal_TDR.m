TimeSec = 3*365*24*60*60;
MACE_Ba_T = 3*1e-04;
WGTS_B_T = 3.6;
BKG_RateSec = 0.01; 
Q_i = 18575;
TD = 'DR30';
range = 30;
Scan = 'ON';
PlotFit = 'OFF'; %only used when Scan ON
SysBudget = '05';
RecomputeFlag = 'ON';
Anchor6G = 'OFF';
%%
LoopArg = {'TimeSec',TimeSec,'Q_i',Q_i,...
    'MACE_Ba_T',MACE_Ba_T,'WGTS_B_T',WGTS_B_T,...
    'TD',TD,'BKG_RateSec',BKG_RateSec,...
    'Scan',Scan,'PlotFit',PlotFit,'SysBudget',SysBudget,...
    'RecomputeFlag',RecomputeFlag,'saveResults','ON','Anchor6G',Anchor6G};
%%
% Do stat + sys (one by one)
 [mNu_Fit, parScan, errScan, chi2minScan, dof ,mNu90,mNumin] = ...
     NuMassScan_SensitivityNominal_Systematics_Loop(LoopArg{:});
 
% Do stat + N+1
[mNu_Fit, parScan, errScan, chi2minScan, dof ,mNu90,mNumin] = ...
    NuMassScan_SensitivityNominal_Nplus1Systematics_Loop(LoopArg{:});

%% Bar Plot
SysBudget = '05';
NuMassScan_Plot_SensitivityNominal_Nplus1Systematics('SingleSys','ON','TD',TD,'SysBudget',SysBudget,...
    'MACE_Ba_T',MACE_Ba_T,'WGTS_B_T',WGTS_B_T,'range',range,'TimeSec',TimeSec,'BKG_RateSec',BKG_RateSec,...
    'SingleSys','ON');
%% RF Breakdown
SysBudget ='05';
NuMassScan_PlotMultiBar_RFBreakdown('TD',TD,'SingleSys','ON','SysBudget',SysBudget,...
   'MACE_Ba_T',MACE_Ba_T,'WGTS_B_T',WGTS_B_T,'range',range,'TimeSec',TimeSec,'BKG_RateSec',BKG_RateSec,...
    'SingleSys','ON');


                  