%% RF Breakdown MultiBar plot
MACE_Ba_T = 3*1e-04;
WGTS_B_T = 1*3.6;
Anchor6G = 'ON';
Anchor6GValue = 3e-3;
Q_i = 18575;
SysBudget = '08';
TimeSec = 365*24*60*60*124/148; 
range = [ 30 , 60];

%% Common Plot Options
TD = 'MTDcreator_E018575.0_30eV_B35_Ba3.0_RedF1.0_NuMF0.40_BkgF0.10_B9';

PlotArg = {'SingleSys','ON',...
   'MACE_Ba_T',MACE_Ba_T,'WGTS_B_T',WGTS_B_T,...
    'SingleSys','ON','TD',TD,'Anchor6G',Anchor6G,'Anchor6GValue',Anchor6GValue};

% MultiBar Plot with main systmatic breakdown
NuMassScan_Plot_SensitivityNominal_Nplus1Systematics(PlotArg{:},'SysBudget',SysBudget,...
    'TimeSec',TimeSec,'TDRSys','OFF','range',range);

%% RF Breakdown MultiBar plot
%NuMassScan_PlotMultiBar_RFBreakdown(PlotArg{:},'SysBudget',SysBudget,'TimeSec',TimeSec,...
%    'range',range);


