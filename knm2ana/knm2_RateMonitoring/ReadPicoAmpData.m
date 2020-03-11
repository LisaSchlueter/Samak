%
% Read Magnus PicoAmpMeter Data
% KNM2 - RW1 Period
% Data Collected by Magnus
% Thierry Lasserre
% Last Modified: 27/02/2020
%

% Folder
foldername = [ './data/'];
dinfo      = dir(fullfile(foldername));
dinfo([dinfo.isdir]) = []; %get rid of all directories

%% Read Current
rwcurr_pathfile = [foldername '/knm2-magnus-RW-current-period1.csv'];
Itmp            = tdfread(rwcurr_pathfile,',');
PicoA_TimeStamp = datetime(Itmp.x0x22Date0x2FTime);  PicoA_TimeStamp.Year = 2019;
PicoA_Current   = Itmp.Rear_Wall_Current_0x5B1100x2DREI0x2D00x2D11100x2D00000x5D0x22*1e9;
PicoA_VoltageEqTmp =  -50 - ( ...
    -0.48634 ...
    + PicoA_Current .* (0.04737) ...
    + PicoA_Current.^2 .* (-0.0013) ...
    + PicoA_Current.^3 .* (9.569e-7))*1e3;
PicoA_VoltageEq = PicoA_VoltageEqTmp - nanmean(PicoA_VoltageEqTmp);

%% Reduce Data by a factor of 10
PicoA_TimeStamp = PicoA_TimeStamp(1:1:end);
PicoA_Current   = PicoA_Current(1:1:end);
PicoA_VoltageEq = PicoA_VoltageEq(1:1:end);


%% Open RW2 Data
RunAnaArg = {'RunList','KNM2_RW1','DataType','Real',...
    'FSDFlag','BlindingKNM2','ELossFlag','KatrinT2',...
    'AnaFlag','StackPixel','RingMerge','Full'};
MR        = MultiRunAnalysis(RunAnaArg{:});
MR.ROIFlag='14keV'; MR.SetROI;
% Time in days
StartTimeStampDays        = days(MR.SingleRunData.StartTimeStamp-MR.SingleRunData.StartTimeStamp(1));
OverallStartTimeStamp     = (MR.SingleRunData.StartTimeStamp);
% Activity
Activity =  (MR.SingleRunData.WGTS_MolFrac_TT+0.5*MR.SingleRunData.WGTS_MolFrac_HT+0.5*MR.SingleRunData.WGTS_MolFrac_DT).*(MR.SingleRunData.WGTS_CD_MolPerCm2);
% HV setting
qUCorrRW   = (MR.SingleRunData.qU_RM(1,:));
qUmeanRW   = mean(MR.SingleRunData.qU_RM(1,:));    
% Include HV Drift
count   = MR.SingleRunData.TBDIS_RM;
sstime  = MR.SingleRunData.qUfrac_RM.*MR.SingleRunData.TimeSec;
rate    = count./sstime;
rateE   = sqrt(count)./sstime;
cf      = MR.RMRateErosCorrectionqUActivity;
rateEquivalentmV     =  -(rate - mean(rate.*cf) ) ./0.7378;
CrateEquivalentmV    =  -(rate.*cf - mean(rate.*cf) ) ./0.7378;
rateEquivalentmV_E   =  (rateE./0.7378);
    

%% Interpolation of RW PicoCurrent/Voltage 
xq = MR.SingleRunData.StartTimeStamp;
PicoA_VoltageEqKNM1 = interp1(PicoA_TimeStamp,PicoA_VoltageEq, xq,'spline');

%% Sanity Plot
figure(2)
subplot(2,1,1)
stairs(xq,CrateEquivalentmV,'LineWidth',3);
hold on
stairs(xq,rateEquivalentmV,'LineWidth',3);
hold off
xlabel('Time');
ylabel('Current (nA)');
PrettyFigureFormat
subplot(2,1,2)
stairs(xq,PicoA_VoltageEqKNM1,'LineWidth',3);
xlim([min(xq) max(xq)])
hold on
stairs(PicoA_TimeStamp,PicoA_VoltageEq,'LineWidth',3);
hold off
xlim([min(xq) max(xq)])
ylim([-50 50])
xlabel('Time');
ylabel('Potential (mV)');
PrettyFigureFormat


%% Summary Plot
figure(5000)
hrw2=stairs(xq,CrateEquivalentmV,'LineWidth',3);
hold on
ha=stairs(xq,PicoA_VoltageEqKNM1,'LineWidth',3);
hold off
xlabel('Time');
ylabel('Potential (mV)');
PrettyFigureFormat
leg=legend([ha hrw2],...
            'Pico-Ampermeter',...
            'RW2 tritium scan data',...
            'location','SouthEast','FontSize',12);
        leg.Color = 'none'; legend boxoff;

