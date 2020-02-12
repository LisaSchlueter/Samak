% Example script how to use the sensitivity class
%% 1. Step: Initialize SensitivityStudy Object (more options available)
TimeSec = (900*24*60*60)*(124/148);
MACE_Ba_T = 7*1e-04;
WGTS_B_T = 0.7*3.6;
Q_i = 18575;
range = 30;
RecomputeFlag = 'OFF';
Scan = 'ON';
PlotFit = 'OFF'; %only used when Scan ON
SysBudget = '08'; % 6=42 day, 07= 300 days, 08 (and manually eloss off)= 900 days
Anchor6G = 'ON';
mNuStop = 1;              % stop value for scan (neutrino mass)
ScanPrcsn = 0.02;
if strcmp(Anchor6G,'OFF')
    FPD_MeanEff = 0.9;%0.95*0.675;
elseif strcmp(Anchor6G,'ON')
    %Anchor6GValue = 335e-3;FPD_MeanEff = 0.95*0.675;
    Anchor6GValue = 500e-3;FPD_MeanEff = 0.95;
end
TD = 'MTDcreator';
chi2 = 'chi2CM';
SysEffect =  'all';
InputArg = {'TimeSec',TimeSec,'Q_i',Q_i,'range',range,...
    'MACE_Ba_T',MACE_Ba_T,'WGTS_B_T',WGTS_B_T,...
    'Scan',Scan,'PlotFit',PlotFit,'SysBudget',SysBudget,'TD',TD,...
    'mNuStop',mNuStop,'ScanPrcsn',ScanPrcsn,...
    'FPD_MeanEff',FPD_MeanEff,'Anchor6G',Anchor6G,'Anchor6GValue',Anchor6GValue,...
    'chi2',chi2,'SysEffect',SysEffect,'RecomputeFlag',RecomputeFlag};

S = SensitivityStudy(InputArg{:});
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
S.PlotSysBreakdownBars('Ranges',30,'SingleSys','ON','TDRSys','OFF','SysInfoBox','ON');
S.PlotRFBreakdownBars('Ranges',[20,30,40,50,70,90],'SingleSys','ON','SysInfoBox','ON');

%% Plot RunTime Evolution
% S.Plot_RunTimeEvolution('TimeSecs',[42 300 600].*124/148*60*60*24,...
%     'SysBudgets',{'06','07','09'},'SingleSys','ON')
%% plot Runtime it for a different range
% S.range = 30;
% S.SetTD; %set TD accordingly
% S.Plot_RunTimeEvolution;

