% Neutrino Mass Sensitivity for first 30 days of KATRIN
TimeDays = 30;
TimeSec = TimeDays*24*60*60*(124/148);
MACE_Ba_T = 7*1e-04;
WGTS_B_T = 0.7*3.6;
MACE_Bmax_T = 0.7*6;
Q_i = 18575;
range = [30, 60];
Scan = 'ON';
PlotFit = 'OFF'; %only used when Scan ON
SysBudget = '06'; % 6=42 day, 07= 300 days, 09 = 900 days
ELossFlag = 'Abdurashitov';
FPD_ROIlow = 14;
FPD_MeanEff = 0.95;
TD = 'MTDcreator';
% 2 Options to give the background
BKG_RateSec = ''; % if none given: computed according to ROI and First Tritium Level
AnchorBkg6G = ''; % if none given: computed according to ROI and First Tritium Level

%% Init Sensitivity Class
Arg = {'TimeSec',TimeSec,...
    'Q_i',Q_i,...
    'range',range(1),...
    'MACE_Ba_T',MACE_Ba_T,...
    'WGTS_B_T',WGTS_B_T,...
    'MACE_Bmax_T',MACE_Bmax_T,...
    'Scan',Scan,...
    'PlotFit',PlotFit,...
    'SysBudget',SysBudget,...
    'TD',TD,...
    'mNuStop',2,...
    'FPD_MeanEff',FPD_MeanEff,...
    'ELossFlag',ELossFlag,...
    'FPD_ROIlow',FPD_ROIlow,...
    'BKG_RateSec',BKG_RateSec,...
    'AnchorBkg6G',AnchorBkg6G,...
    'chi2','chi2Stat'};

S = SensitivityStudy(Arg{:});
 
%% Compute Loops or Look at Data
   S.RecomputeFlag = 'ON'; %in case never computed: switch to ON!
   for i=1:2
       S.range = range(i); S.SetTD; S.InitializeModels;  
       S.NuMassScan_SensitivityNominal_Nplus1Systematics_Loop;
       S.NuMassScan_SystematicsLoop;
   end

%% Plot Breakdowns MultiBar Plots
xLimits = [600, 1000]; % good for 30 days
S.PlotSysBreakdownBars('Ranges', [30,60],'xLimits',xLimits,'RFBreakdown','ON');%,'xLimits',xLimits);
% S.PlotRFBreakdownBars('Ranges', [30,60],'xLimits',xLimits);