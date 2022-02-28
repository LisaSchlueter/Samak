% get knm2 rate monitor points
AnaFlag   = 'StackPixel';

savedir = [getenv('SamakPath'),'knm2ana/knm2_RateMonitoring/results/'];
savename = sprintf('%sknm2_RateMonitoring_GetCounts_%s.mat',savedir,AnaFlag);

if exist(savename,'file') && 1==2
    load(savename);
else
    RunAnaArg = {'DataType','Real',...
        'RadiativeFlag','ON',...
        'minuitOpt','min ; minos',...
        'FSDFlag','KNM2_0p1eV',...
        'ELossFlag','KatrinT2A20',...
        'SysBudget',40,...
        'AnaFlag',AnaFlag,...
        'chi2','chi2Stat',...
        'NonPoissonScaleFactor',1.112,...
        'FSD_Sigma',sqrt(0.0124+0.0025),...
        'RingMerge','None',...
        'PullFlag',99,...;%99 = no pull
        'BKG_PtSlope',3*1e-06,...
        'DopplerEffectFlag','FSD'};
    A = MultiRunAnalysis('RunList','KNM2_Prompt',RunAnaArg{:});
    
    
    %%
    qU               = A.SingleRunData.qU_RM;
    TimeSec_ScanStep = A.SingleRunData.qUfrac_RM.*A.SingleRunData.TimeSec;
    StartTimeStamp = A.SingleRunData.StartTimeStamp;
    LiveTime       = hours(A.SingleRunData.StartTimeStamp-A.SingleRunData.StartTimeStamp(1));
    Counts         = A.SingleRunData.TBDIS_RM;
    Rates          = Counts./ TimeSec_ScanStep;
    RatesErr       = sqrt(Counts)./TimeSec_ScanStep;
    Activity       = A.SingleRunData.WGTS_CD_MolPerCm2.*(A.SingleRunData.WGTS_MolFrac_TT+0.5.*A.SingleRunData.WGTS_MolFrac_DT+0.5.*A.SingleRunData.WGTS_MolFrac_HT);
    
    Idx_Rw1 = abs(A.SingleRunData.RW_BiasVoltage+0.0496)<0.01;
    Idx_Rw2 = abs(A.SingleRunData.RW_BiasVoltage+0.0077)<0.01;
    Idx_Rw3 = abs(A.SingleRunData.RW_BiasVoltage-0.193)<0.01;
    
   
     %Correction_Activity = zeros(1,A.nRuns);
    %     Correction_Activity(Idx_Rw1) = mean(Activity(Idx_Rw1))./Activity(Idx_Rw1);
    %     Correction_Activity(Idx_Rw2) = mean(Activity(Idx_Rw2))./Activity(Idx_Rw2);
    %     Correction_Activity(Idx_Rw3) = mean(Activity(Idx_Rw3))./Activity(Idx_Rw3);
        %     qUSlope_abs1  = dqUCorr.par(1).*(mean(Rates(Idx_Rw1))./mean(dqUCorr.Rate)); % cps/eV
    %     qUSlope_abs2  = dqUCorr.par(1).*(mean(Rates(Idx_Rw2))./mean(dqUCorr.Rate)); % cps/eV
    %     qUSlope_abs3  = dqUCorr.par(1).*(mean(Rates(Idx_Rw3))./mean(dqUCorr.Rate)); % cps/eV
    %     qUDiff1 = mean(qU(Idx_Rw1))-qU(Idx_Rw1);
    %     qUDiff2 = mean(qU(Idx_Rw2))-qU(Idx_Rw2);
    %     qUDiff3 = mean(qU(Idx_Rw3))-qU(Idx_Rw3);
    %     %Correction_HV  = zeros(A.nRuns,1);
    %     Correction_HV(Idx_Rw1) = qUDiff1.*qUSlope_abs1;
    %     Correction_HV(Idx_Rw2) = qUDiff2.*qUSlope_abs2;
%     %     Correction_HV(Idx_Rw3) = qUDiff3.*qUSlope_abs3;
%      RatesCorr  = zeros(A.nRuns,1);
%     RatesCorr(Idx_Rw1) = Rates(Idx_Rw1).*Correction_Activity(Idx_Rw1)+Correction_HV(Idx_Rw1)';
%     RatesCorr(Idx_Rw2) = Rates(Idx_Rw2).*Correction_Activity(Idx_Rw2)+Correction_HV(Idx_Rw2)';
%     RatesCorr(Idx_Rw3) = Rates(Idx_Rw3).*Correction_Activity(Idx_Rw3)+Correction_HV(Idx_Rw3)';
%   
    
    % activity correction
    Correction_Activity = mean(Activity)./Activity;
    
    % retarding potential corrections: load conversion
    dqUCorr = importdata([savedir,'knm2_RateMonitoring_qUSlope.mat']);
    qUSlope_rel = dqUCorr.par(1)./ mean(dqUCorr.Rate); % 1/eV
    qUSlope_abs = dqUCorr.par(1).*(mean(Rates)./mean(dqUCorr.Rate)); % cps/eV
    
    qUDiff = mean(qU)-qU;
    Correction_HV = qUDiff.*qUSlope_abs;
    
    RatesCorr= Rates.*Correction_Activity+Correction_HV;
    RatesCorrErr       = sqrt(RatesCorr.*TimeSec_ScanStep)./TimeSec_ScanStep;
    
    GetFigure
    subplot(3,1,1)
    plot(LiveTime(:),Rates(:))
    hold on;
    plot(LiveTime(:),RatesCorr(:),'Color',rgb('Black'))
    ylabel(sprintf('Rate (cps)'));
    
    subplot(3,1,2)
    plot(LiveTime,(Activity-mean(Activity))./mean(Activity)*1e2);
    ylabel(sprintf('\\Delta Activity (%%)'));
    
    subplot(3,1,3)
    plot(LiveTime,1e3.*qUDiff);
    ylabel(sprintf('\\Delta qU (meV)'));
    % plot(LiveTime,(Rates.*Correction_Activity')+Correction_HV');
    
    %     % convert corrected rates into meV equivalent
    %     %MeVequi = (RatesCorr-mean(RatesCorr))./qUSlope_abs1;
    %
    %     Shifts_Cps = zeros(3,1);
    %     Shifts_Cps(1) = mean(RatesCorr(Idx_Rw1))-mean(RatesCorr);
    %     Shifts_Cps(2) = mean(RatesCorr(Idx_Rw2))-mean(RatesCorr);
    %     Shifts_Cps(3) = mean(RatesCorr(Idx_Rw3))-mean(RatesCorr);
    %
    %     qUSlope_absAll  = dqUCorr.par(1).*mean(RatesCorr)./mean(dqUCorr.Rate); % cps/eV
    %     Shift_meV = 1e3.*Shifts_Cps./qUSlope_absAll;
    %% 
    MakeDir(savedir);
    save(savename,'Rates','RatesCorr','Counts','qU','TimeSec_ScanStep','LiveTime','StartTimeStamp',...
        'Activity','Correction_Activity','Correction_HV',...
        'RatesCorrErr','RatesErr',...
        'Idx_Rw1','Idx_Rw2','Idx_Rw3','dqUCorr')
end
%%

