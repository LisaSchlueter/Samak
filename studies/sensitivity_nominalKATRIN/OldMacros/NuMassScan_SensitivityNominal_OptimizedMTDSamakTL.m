TimeSec = 3*(365*24*60*60)*(124/148);
MACE_Ba_T = 7*1e-04;
WGTS_B_T = 0.7*3.6;
Q_i = 18575;
range = [30];
Scan = 'ON';
PlotFit = 'OFF'; %only used when Scan ON
SysBudget = {'08'}; % 6=42 day, 07= 300 days, 08 (and manually eloss off)= 42 days
RecomputeFlag = 'OFF';
Anchor6G = 'ON';
if strcmp(Anchor6G,'OFF')
    FPD_MeanEff = 0.9;%0.95*0.675;
elseif strcmp(Anchor6G,'ON')
    FPD_MeanEff = 0.95*0.675;
    %FPD_MeanEff = 0.9;
    Anchor6GValue = 335e-3;
end


%%
for s = 1:numel(SysBudget)
    for r=1:numel(range)
        TD  = sprintf('OptimizeMTD%.0f_NuMassFactor_BkgFraction_all_03',range(r));
%        TD  = sprintf('MTDcreator%.0f',range(r));
%TD  = sprintf('MTDcreator_E018575.0_%.0feV_B35_Ba3.0_RedF1.0_NuMF0.40_BkgF0.10_B9',range(r));
        LoopArg = {'TimeSec',TimeSec,'Q_i',Q_i,'range',range(r),...
            'MACE_Ba_T',MACE_Ba_T,'WGTS_B_T',WGTS_B_T,...
            'Scan','OFF','PlotFit',PlotFit,'SysBudget',SysBudget{s},...
            'RecomputeFlag',RecomputeFlag,'saveResults','ON','TD',TD,'mNuStop',1.5,'Anchor6G',Anchor6G,'Anchor6GValue',Anchor6GValue};
        
        % Do stat + sys (one by one)
        NuMassScan_SensitivityNominal_Systematics_Loop(LoopArg{:});        
        % Do stat + N+1
        [mNu_Fit, parScan, errScan, chi2minScan, dof ,mNu90,mNumin] = ...
            NuMassScan_SensitivityNominal_Nplus1Systematics_Loop(LoopArg{:});
    end
end

%% MultiBar Plot with main systmatic breakdown
% SysBudget = '04';
% NuMassScan_Plot_SensitivityNominal_Nplus1Systematics('SingleSys','ON','SysBudget',SysBudget,...
%    'MACE_Ba_T',MACE_Ba_T,'WGTS_B_T',WGTS_B_T,'range',range,'TimeSec',TimeSec,...
%     'SingleSys','ON','TDRSys','OFF','TD','Optimize');
% %% RF Breakdown MultiBar plot
% SysBudget = '03';
% NuMassScan_PlotMultiBar_RFBreakdown('SingleSys','ON','SysBudget',SysBudget,...
%    'MACE_Ba_T',MACE_Ba_T,'WGTS_B_T',WGTS_B_T,'range',range,'TimeSec',TimeSec,...
%     'SingleSys','ON','TD','Optimize');

%% RF Breakdown MultiBar plot
SysBudget = '08';
TimeSec = 3*365*24*60*60*124/148; 
%TimeSec = 42*24*60*60*124/148; 
range = [ 30];
%% Common Plot Options
if strcmp(Anchor6G,'ON')
    TD = 'MTDcreator';
elseif strcmp(Anchor6G,'OFF')
    TD = 'Optimize';
end
PlotArg = {'SingleSys','ON',...
   'MACE_Ba_T',MACE_Ba_T,'WGTS_B_T',WGTS_B_T,...
    'SingleSys','ON','TD',TD,'Anchor6G',Anchor6G,'Anchor6GValue',Anchor6GValue};
% MultiBar Plot with main systmatic breakdown
TimeSec = 3*(365*24*60*60)*(124/148);
SysBudget = '08';
range = [30,60];
NuMassScan_Plot_SensitivityNominal_Nplus1Systematics(PlotArg{:},'SysBudget',SysBudget,...
    'TimeSec',TimeSec,'TDRSys','OFF','range',range);
%% RF Breakdown MultiBar plot
SysBudget = '08';
range = [30,60];
TimeSec = (365*24*60*60)*(124/148);
NuMassScan_PlotMultiBar_RFBreakdown(PlotArg{:},'SysBudget',SysBudget,'TimeSec',TimeSec,...
    'range',range);
%% RunTime Plot
%range = 60;
%NuMassScan_PlotMultiBar_RunTime(PlotArg{:},'SysBudget',{'08_NoELoss','07','06'},...
%     'TimeSec',(124/148*24*60*60).*[900, 300, 42],'range',range);


