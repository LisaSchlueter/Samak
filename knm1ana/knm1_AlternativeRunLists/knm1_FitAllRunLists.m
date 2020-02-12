%
% Perform final Fit with Golden Run List
% and Alternative Run Lists
% Display the results on MC/Data
%
% T. Lasserre
% Last Updated: 12/07/2019
%
clear all;

%% common settings
SysSwitch = 'ON';
exclDataStart         = 14;
nTrials               = 1000;
RecomputeFlag         = 'OFF';
SysEffects            = struct('TASR','ON','FSD','ON','RF_RX','ON','RF_EL','ON','RF_BF','ON','BkgShape','ON','TCoff_RAD','ON','TCoff_OTHER','ON','Stack','ON');
BkgCM                 = 'ON';


%% RunList: KNM1
RunList               = 'KNM1';
MRL = MultiRunAnalysis('RunList',RunList,...
    'chi2','chi2Stat','DataType','Real','exclDataStart',exclDataStart,...
    'fixPar','5 6 7 8 9 10 11','NonPoissonScaleFactor',1.064,...
    'minuitOpt','min ; minos');
% Stat
MRL.chi2='chi2Stat'; MRL.Fit; MRLs = MRL.FitResult;
% CM Shape
switch SysSwitch
    case 'ON'
        MRL.chi2='chi2CMShape';MRL.ComputeCM('SysEffects',SysEffects,'BkgCM',BkgCM);
        MRL.Fit; MRLss = MRL.FitResult;
end

%% RunList: 'KNM1_m149mvRW'
RunList               = 'KNM1_m149mvRW';
MRLm149mvRW = MultiRunAnalysis('RunList',RunList,...
    'chi2','chi2Stat','DataType','Real','exclDataStart',exclDataStart,...
    'fixPar','5 6 7 8 9 10 11','NonPoissonScaleFactor',1.064);
% Stat
MRLm149mvRW.chi2='chi2Stat'; MRLm149mvRW.Fit; MRLm149mvRWs = MRLm149mvRW.FitResult;
% CM Shape
switch SysSwitch
    case 'ON'
        MRLm149mvRW.chi2='chi2CMShape';  MRLm149mvRW.ComputeCM('SysEffects',SysEffects,'BkgCM',BkgCM);
        MRLm149mvRW.Fit; MRLm149mvRWss = MRLm149mvRW.FitResult;
end

%% RunList: 'KNM1-RhoD'
RunList               = 'KNM1-RhoD';
MRLRhoD = MultiRunAnalysis('RunList',RunList,...
    'chi2','chi2Stat','DataType','Real','exclDataStart',exclDataStart,...
    'fixPar','5 6 7 8 9 10 11','NonPoissonScaleFactor',1.064);
% Stat
MRLRhoD.chi2='chi2Stat'; MRLRhoD.Fit; MRLRhoDs = MRLRhoD.FitResult;
% CM Shape
switch SysSwitch
    case 'ON'
        MRLRhoD.chi2='chi2CMShape';MRLRhoD.ComputeCM('SysEffects',SysEffects,'BkgCM',BkgCM);
        MRLRhoD.Fit; MRLRhoDss = MRLRhoD.FitResult;
end

%% RunList: 'KNM1-TpurityLow'
RunList               = 'KNM1-TpurityLow';
MRLTpurityLow = MultiRunAnalysis('RunList',RunList,...
    'chi2','chi2Stat','DataType','Real','exclDataStart',exclDataStart,...
    'fixPar','5 6 7 8 9 10 11','NonPoissonScaleFactor',1.064);
% Stat
MRLTpurityLow.chi2='chi2Stat'; MRLTpurityLow.Fit; MRLTpurityLows = MRLTpurityLow.FitResult;
% CM Shape
switch SysSwitch
    case 'ON'
        MRLTpurityLow.chi2='chi2CMShape';MRLTpurityLow.ComputeCM('SysEffects',SysEffects,'BkgCM',BkgCM);
        MRLTpurityLow.Fit; MRLTpurityLowss = MRLTpurityLow.FitResult;
end

%% RunList: 'KNM1-TpurityHigh'
RunList               = 'KNM1-TpurityHigh';
MRLTpurityHigh = MultiRunAnalysis('RunList',RunList,...
    'chi2','chi2Stat','DataType','Real','exclDataStart',exclDataStart,...
    'fixPar','5 6 7 8 9 10 11','NonPoissonScaleFactor',1.064);
% Stat
MRLTpurityHigh.chi2='chi2Stat'; MRLTpurityHigh.Fit; MRLTpurityHighs = MRLTpurityHigh.FitResult;
% CM Shape
switch SysSwitch
    case 'ON'
        MRLTpurityHigh.chi2='chi2CMShape';MRLTpurityHigh.ComputeCM('SysEffects',SysEffects,'BkgCM',BkgCM);
        MRLTpurityHigh.Fit; MRLTpurityHighss = MRLTpurityHigh.FitResult;
end

%% RunList: 'KNM1-FirstHalfTime'
RunList               = 'KNM1-FirstHalfTime';
MRLFirstHalfTime = MultiRunAnalysis('RunList',RunList,...
    'chi2','chi2Stat','DataType','Real','exclDataStart',exclDataStart,...
    'fixPar','5 6 7 8 9 10 11','NonPoissonScaleFactor',1.064);
% Stat
MRLFirstHalfTime.chi2='chi2Stat'; MRLFirstHalfTime.Fit; MRLFirstHalfTimes = MRLFirstHalfTime.FitResult;
% CM Shape
switch SysSwitch
    case 'ON'
        MRLFirstHalfTime.chi2='chi2CMShape';MRLFirstHalfTime.ComputeCM('SysEffects',SysEffects,'BkgCM',BkgCM);
        MRLFirstHalfTime.Fit; MRLFirstHalfTimess = MRLFirstHalfTime.FitResult;
end

%% RunList: 'KNM1-MiddleHalfTime'
RunList               = 'KNM1-MiddleHalfTime';
MRLMiddleHalfTime = MultiRunAnalysis('RunList',RunList,...
    'chi2','chi2Stat','DataType','Real','exclDataStart',exclDataStart,...
    'fixPar','5 6 7 8 9 10 11','NonPoissonScaleFactor',1.064);
% Stat
MRLMiddleHalfTime.chi2='chi2Stat'; MRLMiddleHalfTime.Fit; MRLMiddleHalfTimes = MRLMiddleHalfTime.FitResult;
% CM Shape
switch SysSwitch
    case 'ON'
        MRLMiddleHalfTime.chi2='chi2CMShape';MRLMiddleHalfTime.ComputeCM('SysEffects',SysEffects,'BkgCM',BkgCM);
        MRLMiddleHalfTime.Fit; MRLMiddleHalfTimess = MRLMiddleHalfTime.FitResult;
end

%% RunList: 'KNM1-LastHalfTime'
RunList               = 'KNM1-LastHalfTime';
MRLLastHalfTime = MultiRunAnalysis('RunList',RunList,...
    'chi2','chi2Stat','DataType','Real','exclDataStart',exclDataStart,...
    'fixPar','5 6 7 8 9 10 11','NonPoissonScaleFactor',1.064);
% Stat
MRLLastHalfTime.chi2='chi2Stat'; MRLLastHalfTime.Fit; MRLLastHalfTimes = MRLLastHalfTime.FitResult;
% CM Shape
switch SysSwitch
    case 'ON'
        MRLLastHalfTime.chi2='chi2CMShape';MRLLastHalfTime.ComputeCM('SysEffects',SysEffects,'BkgCM',BkgCM);
        MRLLastHalfTime.Fit; MRLLastHalfTimess = MRLLastHalfTime.FitResult;
end



%% RunList: 'KNM1downScan'
RunList               = 'KNM1downScan';
MRLdownScan = MultiRunAnalysis('RunList',RunList,...
    'chi2','chi2Stat','DataType','Real','exclDataStart',exclDataStart,...
    'fixPar','5 6 7 8 9 10 11','NonPoissonScaleFactor',1.064);
% Stat
MRLdownScan.chi2='chi2Stat'; MRLdownScan.Fit; MRLdownScans = MRLdownScan.FitResult;
% CM Shape
switch SysSwitch
    case 'ON'
        MRLdownScan.chi2='chi2CMShape';MRLdownScan.ComputeCM('SysEffects',SysEffects,'BkgCM',BkgCM);
        MRLdownScan.Fit; MRLdownScanss = MRLdownScan.FitResult;
end



%% RunList: 'KNM1upScan'
RunList               = 'KNM1upScan';
MRLupScan = MultiRunAnalysis('RunList',RunList,...
    'chi2','chi2Stat','DataType','Real','exclDataStart',exclDataStart,...
    'fixPar','5 6 7 8 9 10 11','NonPoissonScaleFactor',1.064);
% Stat
MRLupScan.chi2='chi2Stat'; MRLupScan.Fit; MRLupScans = MRLupScan.FitResult;
% CM Shape
switch SysSwitch
    case 'ON'
        MRLupScan.chi2='chi2CMShape';MRLupScan.ComputeCM('SysEffects',SysEffects,'BkgCM',BkgCM);
        MRLupScan.Fit; MRLupScanss = MRLupScan.FitResult;
end



%% Display Results - Stat Only
t = PrintTable('Sensitivity Fit Results with Various RunLists');
t.addRow('Run List','m^2_\beta (eV^2)','$E_0$ (eV)','B (mcps)','$\chi^2$ / dof');

t.addRow(sprintf('%s ',MRL.StackFileName),...
    sprintf('%.3f \\pm %.3f',MRLs.par(1), MRLs.err(1)),...
    sprintf('%.2f \\pm %.2f',MRL.ModelObj.Q_i+MRLs.par(2),MRLs.err(2)),...
    sprintf('%.1f \\pm %.1f',(MRL.ModelObj.BKG_RateSec_i + MRLs.par(3))*1e3,MRLs.err(3)*1e3),...
    sprintf('%.1f / %.0f dof',MRLs.chi2min,MRLs.dof));

t.addRow(sprintf('%s ',MRLm149mvRW.StackFileName),...
    sprintf('%.3f \\pm %.3f',MRLm149mvRWs.par(1), MRLm149mvRWs.err(1)),...
    sprintf('%.2f \\pm %.2f',MRLm149mvRW.ModelObj.Q_i+MRLm149mvRWs.par(2),MRLm149mvRWs.err(2)),...
    sprintf('%.1f \\pm %.1f',(MRLm149mvRW.ModelObj.BKG_RateSec_i + MRLm149mvRWs.par(3))*1e3,MRLm149mvRWs.err(3)*1e3),...
    sprintf('%.1f / %.0f dof',MRLm149mvRWs.chi2min,MRLm149mvRWs.dof));

t.addRow(sprintf('%s ',MRLRhoD.StackFileName),...
    sprintf('%.3f \\pm %.3f',MRLRhoDs.par(1), MRLRhoDs.err(1)),...
    sprintf('%.2f \\pm %.2f',MRLRhoD.ModelObj.Q_i+MRLRhoDs.par(2),MRLRhoDs.err(2)),...
    sprintf('%.1f \\pm %.1f',(MRLRhoD.ModelObj.BKG_RateSec_i + MRLRhoDs.par(3))*1e3,MRLRhoDs.err(3)*1e3),...
    sprintf('%.1f / %.0f dof',MRLRhoDs.chi2min,MRLRhoDs.dof));

t.addRow(sprintf('%s ',MRLTpurityLow.StackFileName),...
    sprintf('%.3f \\pm %.3f',MRLTpurityLows.par(1), MRLTpurityLows.err(1)),...
    sprintf('%.2f \\pm %.2f',MRLTpurityLow.ModelObj.Q_i+MRLTpurityLows.par(2),MRLTpurityLows.err(2)),...
    sprintf('%.1f \\pm %.1f',(MRLTpurityLow.ModelObj.BKG_RateSec_i + MRLTpurityLows.par(3))*1e3,MRLTpurityLows.err(3)*1e3),...
    sprintf('%.1f / %.0f dof',MRLTpurityLows.chi2min,MRLTpurityLows.dof));

t.addRow(sprintf('%s ',MRLTpurityHigh.StackFileName),...
    sprintf('%.3f \\pm %.3f',MRLTpurityHighs.par(1), MRLTpurityHighs.err(1)),...
    sprintf('%.2f \\pm %.2f',MRLTpurityHigh.ModelObj.Q_i+MRLTpurityHighs.par(2),MRLTpurityHighs.err(2)),...
    sprintf('%.1f \\pm %.1f',(MRLTpurityHigh.ModelObj.BKG_RateSec_i + MRLTpurityHighs.par(3))*1e3,MRLTpurityHighs.err(3)*1e3),...
    sprintf('%.1f / %.0f dof',MRLTpurityHighs.chi2min,MRLTpurityHighs.dof));

t.addRow(sprintf('%s ',MRLFirstHalfTime.StackFileName),...
    sprintf('%.3f \\pm %.3f',MRLFirstHalfTimes.par(1), MRLFirstHalfTimes.err(1)),...
    sprintf('%.2f \\pm %.2f',MRLFirstHalfTime.ModelObj.Q_i+MRLFirstHalfTimes.par(2),MRLFirstHalfTimes.err(2)),...
    sprintf('%.1f \\pm %.1f',(MRLFirstHalfTime.ModelObj.BKG_RateSec_i + MRLFirstHalfTimes.par(3))*1e3,MRLFirstHalfTimes.err(3)*1e3),...
    sprintf('%.1f / %.0f dof',MRLFirstHalfTimes.chi2min,MRLFirstHalfTimes.dof));

t.addRow(sprintf('%s ',MRLMiddleHalfTime.StackFileName),...
    sprintf('%.3f \\pm %.3f',MRLMiddleHalfTimes.par(1), MRLMiddleHalfTimes.err(1)),...
    sprintf('%.2f \\pm %.2f',MRLMiddleHalfTime.ModelObj.Q_i+MRLMiddleHalfTimes.par(2),MRLMiddleHalfTimes.err(2)),...
    sprintf('%.1f \\pm %.1f',(MRLMiddleHalfTime.ModelObj.BKG_RateSec_i + MRLMiddleHalfTimes.par(3))*1e3,MRLMiddleHalfTimes.err(3)*1e3),...
    sprintf('%.1f / %.0f dof',MRLMiddleHalfTimes.chi2min,MRLMiddleHalfTimes.dof));

t.addRow(sprintf('%s ',MRLLastHalfTime.StackFileName),...
    sprintf('%.3f \\pm %.3f',MRLLastHalfTimes.par(1), MRLLastHalfTimes.err(1)),...
    sprintf('%.2f \\pm %.2f',MRLLastHalfTime.ModelObj.Q_i+MRLLastHalfTimes.par(2),MRLLastHalfTimes.err(2)),...
    sprintf('%.1f \\pm %.1f',(MRLLastHalfTime.ModelObj.BKG_RateSec_i + MRLLastHalfTimes.par(3))*1e3,MRLLastHalfTimes.err(3)*1e3),...
    sprintf('%.1f / %.0f dof',MRLLastHalfTimes.chi2min,MRLLastHalfTimes.dof));

t.addRow(sprintf('%s ',MRLdownScan.StackFileName),...
    sprintf('%.3f \\pm %.3f',MRLdownScans.par(1), MRLdownScans.err(1)),...
    sprintf('%.2f \\pm %.2f',MRLdownScan.ModelObj.Q_i+MRLdownScans.par(2),MRLdownScans.err(2)),...
    sprintf('%.1f \\pm %.1f',(MRLdownScan.ModelObj.BKG_RateSec_i + MRLdownScans.par(3))*1e3,MRLdownScans.err(3)*1e3),...
    sprintf('%.1f / %.0f dof',MRLdownScans.chi2min,MRLdownScans.dof));

t.addRow(sprintf('%s ',MRLupScan.StackFileName),...
    sprintf('%.3f \\pm %.3f',MRLupScans.par(1), MRLupScans.err(1)),...
    sprintf('%.2f \\pm %.2f',MRLupScan.ModelObj.Q_i+MRLupScans.par(2),MRLupScans.err(2)),...
    sprintf('%.1f \\pm %.1f',(MRLupScan.ModelObj.BKG_RateSec_i + MRLupScans.par(3))*1e3,MRLupScans.err(3)*1e3),...
    sprintf('%.1f / %.0f dof',MRLupScans.chi2min,MRLupScans.dof));

t.display;
t.HasHeader = true;
t.Format = 'tex';
t.Caption = sprintf('Fit Results  - Fitter = %s - Range = [%.1f - %.1f] eV - Statistical Error Only',...
    MRL.fitter,MRL.ModelObj.qU(MRL.exclDataStart),MRL.ModelObj.qU(end));
t.print

%% Display Results Stat+Syst
t = PrintTable('Sensitivity Fit Results with Various RunLists');
t.addRow('Run List','m^2_\beta (eV^2)','$E_0$ (eV)','B (mcps)','$\chi^2$ / dof');

t.addRow(sprintf('%s ',MRL.StackFileName),...
    sprintf('%.3f \\pm %.3f',MRLss.par(1), MRLss.err(1)),...
    sprintf('%.2f \\pm %.2f',MRL.ModelObj.Q_i+MRLss.par(2),MRLss.err(2)),...
    sprintf('%.1f \\pm %.1f',(MRL.ModelObj.BKG_RateSec_i + MRLss.par(3))*1e3,MRLss.err(3)*1e3),...
    sprintf('%.1f / %.0f dof',MRLss.chi2min,MRLss.dof));

t.addRow(sprintf('%s ',MRLm149mvRW.StackFileName),...
    sprintf('%.3f \\pm %.3f',MRLm149mvRWss.par(1), MRLm149mvRWss.err(1)),...
    sprintf('%.2f \\pm %.2f',MRLm149mvRW.ModelObj.Q_i+MRLm149mvRWss.par(2),MRLm149mvRWss.err(2)),...
    sprintf('%.1f \\pm %.1f',(MRLm149mvRW.ModelObj.BKG_RateSec_i + MRLm149mvRWss.par(3))*1e3,MRLm149mvRWss.err(3)*1e3),...
    sprintf('%.1f / %.0f dof',MRLm149mvRWss.chi2min,MRLm149mvRWss.dof));

t.addRow(sprintf('%s ',MRLRhoD.StackFileName),...
    sprintf('%.3f \\pm %.3f',MRLRhoDss.par(1), MRLRhoDss.err(1)),...
    sprintf('%.2f \\pm %.2f',MRLRhoD.ModelObj.Q_i+MRLRhoDss.par(2),MRLRhoDss.err(2)),...
    sprintf('%.1f \\pm %.1f',(MRLRhoD.ModelObj.BKG_RateSec_i + MRLRhoDss.par(3))*1e3,MRLRhoDss.err(3)*1e3),...
    sprintf('%.1f / %.0f dof',MRLRhoDss.chi2min,MRLRhoDss.dof));

t.addRow(sprintf('%s ',MRLTpurityLow.StackFileName),...
    sprintf('%.3f \\pm %.3f',MRLTpurityLowss.par(1), MRLTpurityLowss.err(1)),...
    sprintf('%.2f \\pm %.2f',MRLTpurityLow.ModelObj.Q_i+MRLTpurityLowss.par(2),MRLTpurityLowss.err(2)),...
    sprintf('%.1f \\pm %.1f',(MRLTpurityLow.ModelObj.BKG_RateSec_i + MRLTpurityLowss.par(3))*1e3,MRLTpurityLowss.err(3)*1e3),...
    sprintf('%.1f / %.0f dof',MRLTpurityLowss.chi2min,MRLTpurityLowss.dof));

t.addRow(sprintf('%s ',MRLTpurityHigh.StackFileName),...
    sprintf('%.3f \\pm %.3f',MRLTpurityHighss.par(1), MRLTpurityHighss.err(1)),...
    sprintf('%.2f \\pm %.2f',MRLTpurityHigh.ModelObj.Q_i+MRLTpurityHighss.par(2),MRLTpurityHighss.err(2)),...
    sprintf('%.1f \\pm %.1f',(MRLTpurityHigh.ModelObj.BKG_RateSec_i + MRLTpurityHighss.par(3))*1e3,MRLTpurityHighss.err(3)*1e3),...
    sprintf('%.1f / %.0f dof',MRLTpurityHighss.chi2min,MRLTpurityHighss.dof));

t.addRow(sprintf('%s ',MRLFirstHalfTime.StackFileName),...
    sprintf('%.3f \\pm %.3f',MRLFirstHalfTimess.par(1), MRLFirstHalfTimess.err(1)),...
    sprintf('%.2f \\pm %.2f',MRLFirstHalfTime.ModelObj.Q_i+MRLFirstHalfTimess.par(2),MRLFirstHalfTimess.err(2)),...
    sprintf('%.1f \\pm %.1f',(MRLFirstHalfTime.ModelObj.BKG_RateSec_i + MRLFirstHalfTimess.par(3))*1e3,MRLFirstHalfTimess.err(3)*1e3),...
    sprintf('%.1f / %.0f dof',MRLFirstHalfTimess.chi2min,MRLFirstHalfTimess.dof));

t.addRow(sprintf('%s ',MRLMiddleHalfTime.StackFileName),...
    sprintf('%.3f \\pm %.3f',MRLMiddleHalfTimess.par(1), MRLMiddleHalfTimess.err(1)),...
    sprintf('%.2f \\pm %.2f',MRLMiddleHalfTime.ModelObj.Q_i+MRLMiddleHalfTimess.par(2),MRLMiddleHalfTimess.err(2)),...
    sprintf('%.1f \\pm %.1f',(MRLMiddleHalfTime.ModelObj.BKG_RateSec_i + MRLMiddleHalfTimess.par(3))*1e3,MRLMiddleHalfTimess.err(3)*1e3),...
    sprintf('%.1f / %.0f dof',MRLMiddleHalfTimess.chi2min,MRLMiddleHalfTimess.dof));

t.addRow(sprintf('%s ',MRLLastHalfTime.StackFileName),...
    sprintf('%.3f \\pm %.3f',MRLLastHalfTimess.par(1), MRLLastHalfTimess.err(1)),...
    sprintf('%.2f \\pm %.2f',MRLLastHalfTime.ModelObj.Q_i+MRLLastHalfTimess.par(2),MRLLastHalfTimess.err(2)),...
    sprintf('%.1f \\pm %.1f',(MRLLastHalfTime.ModelObj.BKG_RateSec_i + MRLLastHalfTimess.par(3))*1e3,MRLLastHalfTimess.err(3)*1e3),...
    sprintf('%.1f / %.0f dof',MRLLastHalfTimess.chi2min,MRLLastHalfTimess.dof));

t.addRow(sprintf('%s ',MRLdownScan.StackFileName),...
    sprintf('%.3f \\pm %.3f',MRLdownScanss.par(1), MRLdownScanss.err(1)),...
    sprintf('%.2f \\pm %.2f',MRLdownScan.ModelObj.Q_i+MRLdownScanss.par(2),MRLdownScanss.err(2)),...
    sprintf('%.1f \\pm %.1f',(MRLdownScan.ModelObj.BKG_RateSec_i + MRLdownScanss.par(3))*1e3,MRLdownScanss.err(3)*1e3),...
    sprintf('%.1f / %.0f dof',MRLdownScanss.chi2min,MRLdownScanss.dof));

t.addRow(sprintf('%s ',MRLupScan.StackFileName),...
    sprintf('%.3f \\pm %.3f',MRLupScanss.par(1), MRLupScanss.err(1)),...
    sprintf('%.2f \\pm %.2f',MRLupScan.ModelObj.Q_i+MRLupScanss.par(2),MRLupScanss.err(2)),...
    sprintf('%.1f \\pm %.1f',(MRLupScan.ModelObj.BKG_RateSec_i + MRLupScanss.par(3))*1e3,MRLupScanss.err(3)*1e3),...
    sprintf('%.1f / %.0f dof',MRLupScanss.chi2min,MRLupScanss.dof));

t.display;
t.HasHeader = true;
t.Format = 'tex';
t.Caption = sprintf('Fit Results  - Fitter = %s - Range = [%.1f - %.1f] eV - Statistical Error Only',...
    MRL.fitter,MRL.ModelObj.qU(MRL.exclDataStart),MRL.ModelObj.qU(end));
t.print
