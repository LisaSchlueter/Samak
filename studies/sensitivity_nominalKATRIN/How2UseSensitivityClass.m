% Example script how to use the sensitivity class
%% 1. Step: Settings
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

%% 2. Step: Init Sensitivity Class
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

%% Example 1: Compute Sensitivity for your settings and Plot
S.NuMassScan; 
S.PlotScan;
%% Example 2: Loop over systematics and display results
S.NuMassScan_SystematicsLoop; % Single Contributions
S.MultimNu90 % Sensitivity on Neutrino Mass Squared 90% C.L. 
%% Example 3: Loop over systematics (stacked) and display results
S.NuMassScan_SensitivityNominal_Nplus1Systematics_Loop; % Stacked: stat, stat+tc, stat+tc+fsd, stat+tc+fsd+rf(all)
S.MultimNu90Stack % Sensitivity on Neutrino Mass Squared 90% C.L. 
%% Plot options: Systematic Breakdowns
% Caution: If result files are not available: they are computed!
S.PlotSysBreakdownBars('Ranges',[30, 60],'SingleSys','ON','TDRSys','ON','SysInfoBox','ON');
S.PlotRFBreakdownBars('Ranges',[30, 60],'SingleSys','ON','SysInfoBox','ON');

%% Plot RunTime Evolution
S.Plot_RunTimeEvolution('TimeSecs',[42 300 600].*124/148*60*60*24,...
    'SysBudgets',{'06','07','09'},'SingleSys','ON')
%% plot Runtime it for a different range
S.range = 60;
S.SetTD; %set TD accordingly
S.Plot_RunTimeEvolution;

