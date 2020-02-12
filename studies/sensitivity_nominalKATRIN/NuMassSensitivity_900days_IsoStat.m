% Neutrino Mass Sensitivity for first 30 days of KATRIN
TimeDays = 900;
TimeSec = TimeDays*24*60*60*(124/148);
MACE_Ba_T = 7*1e-04;
WGTS_B_T = 0.7*3.6;
MACE_Bmax_T = 0.7*6;
Q_i = 18575;
range = [30,60];
Scan = 'ON';
PlotFit = 'OFF'; %only used when Scan ON
SysBudget = '09'; % 6=42 day, 07= 300 days, 09 = 900 days
ELossFlag = 'Abdurashitov';
FPD_ROIlow = 14;
FPD_MeanEff = 0.95;
BKG_RateSec = ''; % if none given: computed according to ROI
TD = 'MTDcreator';
%AnchorBkg6G = '';
AnchorBkg6G = 0.1;

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
    'mNuStop',1,...
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
%xLimits = 
S.PlotSysBreakdownBars('Ranges', [30,60],'xLimits',xLimits,'RFBreakdown','ON');%,'xLimits',xLimits);
%S.PlotSysBreakdownBars('Ranges', [90,120],'xLimits',xLimits,'RFBreakdown','ON');%,'xLimits',xLimits);
[280, 350]; % good for 900 days
xLimits = [250, 340]; % good for 900 daysS.PlotRFBreakdownBars('Ranges', [30,60],'xLimits',xLimits);