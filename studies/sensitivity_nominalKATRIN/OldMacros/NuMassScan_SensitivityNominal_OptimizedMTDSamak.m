TimeSec = 30*24*60*60*(124/148);%(3*365*24*60*60);%(36*24*60*60);%*(124/148);
MACE_Ba_T = 7*1e-04;
WGTS_B_T = 0.7*3.6;
MACE_Bmax_T = 0.7*6;
Q_i = 18575;
range = [30, 60];
Scan = 'ON';
PlotFit = 'OFF'; %only used when Scan ON
SysBudget = '06'; % 6=42 day, 07= 300 days, 09 = 900 days
RecomputeFlag = 'OFF';
ELossFlag = 'Abdurashitov';
FPD_ROIlow = 14;
FPD_MeanEff = 0.95;%5;
BKG_RateSec = ''; % if none given: computed according to ROI
TD = 'MTDcreator';
AnchorBkg6G = '';%[0.05907,0.1182,0.2363];
%for b = 1:numel(AnchorBkg6G_all)
%AnchorBkg6G = AnchorBkg6G_all(b);
%% Common Plot Options

Arg = {'TimeSec',TimeSec,'Q_i',Q_i,'range',range(1),...
    'MACE_Ba_T',MACE_Ba_T,'WGTS_B_T',WGTS_B_T,'MACE_Bmax_T',MACE_Bmax_T,...
    'Scan',Scan,'PlotFit',PlotFit,'SysBudget',SysBudget,...
    'RecomputeFlag',RecomputeFlag,'TD',TD,'mNuStop',2,...
    'FPD_MeanEff',FPD_MeanEff,...
    'ELossFlag',ELossFlag,'FPD_ROIlow',FPD_ROIlow,...
    'BKG_RateSec',BKG_RateSec,'AnchorBkg6G',AnchorBkg6G};

S = SensitivityStudy(Arg{:});
 
%% Compute Loops or Look at Data
% % %S.range = 60; S.SetTD;
   S.RecomputeFlag = 'OFF';
   for i=1:2
       S.range = range(i); S.SetTD; S.InitializeModels;  
       S.NuMassScan_SensitivityNominal_Nplus1Systematics_Loop;
       S.NuMassScan_SystematicsLoop;
   end

% %  
% fprintf('-----Neurino Mass sensitivities:------ \n')
% fprintf('Stat = %.0f meV\n +TC = %.0f meV \n+FSD = %.0f meV \n +RF = %.0f meV \n',...
%     sqrt(struct2array(S.MultimNu90Stack))*1e3)
%% Plot Breakdowns MultiBar Plots
%42 days ---'
switch TimeSec./(24*60*60.*124/148)
    case 30
        xLimits = [680, 820]; % 30 day
    case 42
        xLimits = [620, 860]; % 42 day
    case 45
        xLimits = [610, 750]; % 42 day
    case 300
        xLimits = [380 535];
    case 900
        xLimits = [280, 350];
end
  xLimits = [680, 1000]; % 30 day
%xLimits = [200, 295]; 
%xLimits = [380 535];
%xLimits = '';
S.PlotSysBreakdownBars('Ranges', [30,60],'xLimits',xLimits);%,'xLimits',xLimits);
S.PlotRFBreakdownBars('Ranges', [30,60],'xLimits',xLimits);
%end
%% RunTime Plot
%S.range = 30; S.SetTD;
%S.Plot_RunTimeEvolution('SysBudgets',{'09_NoELoss','07','06'},'TimeSec',(124/148*24*60*60).*[900, 300, 45]);
%S.Plot_RunTimeEvolution('SysBudgets',{'09_NoELoss','07','06','06','06'},'TimeSec',(124/148*24*60*60).*[900,300,60, 45, 30]);