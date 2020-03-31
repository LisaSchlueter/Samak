KNM1CorFlag    = 'ON';
HVdriftCorFlag = 'OFF';
ROI            = '14keV';

%% KNM1 with Calibration - Divide Rates by
switch KNM1CorFlag
    case 'OFF'
        KNM1correction  = [ 1.0000    1          1         1     ];
    case 'ON'
        KNM1correction  = [ 1.0000    0.9992    0.9975    0.9961];
end

% HV Drift
switch HVdriftCorFlag
    case 'ON'
        HVdriftMvDay = 1.5;
    case 'OFF'
        HVdriftMvDay = 0.0;
end

%
% Reference Rate: KNM2_RW2
%
    DataType  = 'Real';
    RunAnaArg = {'RunList','KNM2_RW2','DataType',DataType,...
        'FSDFlag','BlindingKNM2','ELossFlag','KatrinT2',...
        'AnaFlag','StackPixel','RingMerge','Full'};
    MR        = MultiRunAnalysis(RunAnaArg{:});
    A         = RingAnalysis('RunAnaObj',MR,'RingList',1:4);
    R         = A.MultiObj(1);
    R.ROIFlag=ROI; R.SetROI;
    

    %% Time in days
    StartTimeStampDays           = days(R.SingleRunData.StartTimeStamp-R.SingleRunData.StartTimeStamp(1));
    
    %% HV Drift Correction - MOS
    FirstDayPeriod1 = datetime('02-Oct-2019 14:52:19');
    TimeLineDaysFirstDayPeriod1 = days(MR.SingleRunData.StartTimeStamp-FirstDayPeriod1);
    HVdriftPerPixel = HVdriftMvDay*(TimeLineDaysFirstDayPeriod1) * 6.303e-3;

    %% Stacked Pixel Data for each patch
    count  = zeros(A.nRings,numel(A.RunAnaObj.RunList));
    rate   = zeros(A.nRings,numel(A.RunAnaObj.RunList));
    rateE  = zeros(A.nRings,numel(A.RunAnaObj.RunList));
    sstime = zeros(A.nRings,numel(A.RunAnaObj.RunList));
    RefPeriodRW2CPS = zeros(A.nRings,1);
    qUmeanRW2       = zeros(A.nRings,1);
    
    for i=1:A.nRings
        R           = A.MultiObj(i);
        R.ROIFlag=ROI; R.SetROI;
        
        % Slow Control Data: HV setting
        qUCorrRW2       = (R.SingleRunData.qU_RM(1,:));
        qUmeanRW2(i)    = mean(qUCorrRW2);
        
        % Correct for activity and qU
        R.RMCorrection('saveplot','OFF','pixlist',sprintf('ring%i',i),'QAplots','OFF');
    
        % Include HV Drift if any
        count(i,:)  = R.SingleRunData.TBDIS_RM + HVdriftPerPixel.*numel(R.PixList).*R.SingleRunData.qUfrac_RM.*R.SingleRunData.TimeSec;
        
        % Correct for KNM1 Radial Effect if any
        count(i,:)  = count(i,:) ./ KNM1correction(i);
        
        sstime(i,:) = R.SingleRunData.qUfrac_RM.*R.SingleRunData.TimeSec;
        rate(i,:)   = count(i,:)./sstime(i,:);
        rateE(i,:)  = sqrt(count(i,:))./sstime(i,:);
        RefPeriodRW2CPS(i) = mean(rate(i,:));
    end
    ActivityRW2  =  mean(R.SingleRunData.WGTS_MolFrac_TT'+0.5*R.SingleRunData.WGTS_MolFrac_HT'+0.5*R.SingleRunData.WGTS_MolFrac_DT')...
        .*mean(R.SingleRunData.WGTS_CD_MolPerCm2);

    %
    % End of Reference Rate: KNM2_RW2
    %
    
    %
    % Loop on Period RW1 RW2 RW3
    %
    
    %%
for j=1:3
    
    RunList   = ['KNM2_RW' num2str(j)];

    % Read Data
    DataType  = 'Real';
    RunAnaArg = {'RunList',RunList,'DataType',DataType,...
        'FSDFlag','BlindingKNM2','ELossFlag','KatrinT2',...
        'AnaFlag','StackPixel','RingMerge','Full'};
    MR        = MultiRunAnalysis(RunAnaArg{:});
    A         = RingAnalysis('RunAnaObj',MR,'RingList',1:4);
    R         = A.MultiObj(1);
    R.ROIFlag=ROI; R.SetROI;

    % Time in days
    StartTimeStampDays           = days(R.SingleRunData.StartTimeStamp-R.SingleRunData.StartTimeStamp(1));
    OverallStartTimeStamp{j}     = (R.SingleRunData.StartTimeStamp);
    
    % HV Drift Correction
    FirstDayPeriod1 = datetime('02-Oct-2019 14:52:19');
    TimeLineDaysFirstDayPeriod1 = days(MR.SingleRunData.StartTimeStamp-FirstDayPeriod1);
    HVdriftPerPixel  = HVdriftMvDay*(TimeLineDaysFirstDayPeriod1) * 6.303e-3;
    
    % Activity
    Activity{j} =  (R.SingleRunData.WGTS_MolFrac_TT+0.5*R.SingleRunData.WGTS_MolFrac_HT+0.5*R.SingleRunData.WGTS_MolFrac_DT)...
        .*(R.SingleRunData.WGTS_CD_MolPerCm2);
    
    % Stacked Pixel Data for each patch
    count   = zeros(A.nRings,numel(A.RunAnaObj.RunList));
    rateRaw = zeros(A.nRings,numel(A.RunAnaObj.RunList));
    rate    = zeros(A.nRings,numel(A.RunAnaObj.RunList));
    rateE   = zeros(A.nRings,numel(A.RunAnaObj.RunList));
    cf      = zeros(A.nRings,numel(A.RunAnaObj.RunList));
    sstime  = zeros(A.nRings,numel(A.RunAnaObj.RunList));
    
    for i=1:A.nRings
        R           = A.MultiObj(i);
        R.ROIFlag   = ROI; 
        R.SetROI;

        % HV setting
        qUCorrRW        = (R.SingleRunData.qU_RM(1,:));
        qUmeanRW{j,i}   = mean(R.SingleRunData.qU_RM(1,:));
        
        sstime(i,:) = R.SingleRunData.qUfrac_RM.*R.SingleRunData.TimeSec;
    
        % Include HV Drift
        count(i,:)  = R.SingleRunData.TBDIS_RM + HVdriftPerPixel.*numel(R.PixList).*R.SingleRunData.qUfrac_RM.*R.SingleRunData.TimeSec;
        
        % Correct for KNM1 Radial Effect
        count(i,:)  = count(i,:)./KNM1correction(i);
        
        % Rate uncorrected for activity and qU
        rateRaw(i,:) = count(i,:)./sstime(i,:);
        
        % Ref = Period2 - Average qU correction  
        HVcorrCPSperPixel{j,i} = (qUmeanRW{j,i} - qUmeanRW2(i)) * 6.3032 * numel(R.PixList);
        
        % Correct for activity and qU 1
        cf(i,:)     = R.RMRateErosCorrectionqUActivity;
        Crate1{j,i} = (rateRaw(i,:).*cf(i,:)./mean(Activity{j}).*ActivityRW2 + HVcorrCPSperPixel{j,i} - RefPeriodRW2CPS(1)./numel(A.MultiObj(1).PixList).*numel(R.PixList));
        
        % Correct for activity and qU
        R.RMCorrection('saveplot','OFF','pixlist',sprintf('ring%i',i),'QAplots','OFF');
        count(i,:)  = R.SingleRunData.TBDIS_RM + HVdriftPerPixel.*numel(R.PixList).*R.SingleRunData.qUfrac_RM.*R.SingleRunData.TimeSec;
        
        rate(i,:)   = count(i,:)./sstime(i,:);
        rateE(i,:)  = sqrt(count(i,:))./sstime(i,:);
        
        
        
        % Ref = Period2 PSR1
        % j = RW Period
        % i = PSR
        Crate2{j,i}                =  (rate(i,:)./mean(Activity{j}).*ActivityRW2 + HVcorrCPSperPixel{j,i} - RefPeriodRW2CPS(1)./numel(A.MultiObj(1).PixList).*numel(R.PixList));
        
        CrateEquivalentmV{j,i}    =  -(rate(i,:)./mean(Activity{j}).*ActivityRW2 + HVcorrCPSperPixel{j,i} - RefPeriodRW2CPS(1)./numel(A.MultiObj(1).PixList).*numel(R.PixList)) ./737.8 * 1e3 * 117 / numel(R.PixList);
        
        rateEquivalentmV_E{j,i}   =  (rateE(i,:) ./737.8 *1e3 * 117 / numel(R.PixList));
        
    end
    
    %% Rate Evolution --> mV equivalent
    myMainTitle = sprintf('KATRIN - KNM2 RW%s - FPD Rate Evolution 300eV below Endpoint',num2str(j));
    maintitle   = myMainTitle;
    savefile1    = sprintf('plots/KNM2_RM300_EffectivePotentialFluctuation_RW%.0f_PseudoRings_1.png',j);
    fig1      = figure('Name',sprintf('KATRIN - %s Scanwise Background',RunList),...
        'NumberTitle','off','rend','painters','pos',[10 10 1400 1000]);
    a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
    a.FontSize=24;a.FontWeight='bold';
    
    for i=1:A.nRings

        subplot(A.nRings,1,i)
        hi=plot(StartTimeStampDays,CrateEquivalentmV{j,i},'s--','Color',rgb('IndianRed'),'LineWidth',1,'MarkerSize',12,'markerfacecolor',rgb('IndianRed'));
        
        % Fit: Constant
        fitType = fittype('a*x+b'); % The equation for your fit goes here
        [f,gofCons] = fit(StartTimeStampDays',CrateEquivalentmV{j,i}',fitType,...
            'Weights',(1./rateEquivalentmV_E{j,i}).^2,...
            'StartPoint',[0 CrateEquivalentmV{j,i}(1)],...
            'Robust','Bisquare',...
            'Lower',[0 -1e9],...
            'Upper',[0 1e9]);
        ci = confint(f,0.68); uncertaintyCons = (ci(2,1)-ci(1,1))/2;
        hold on
        hc = plot(StartTimeStampDays,f.a.*ones(1,numel(StartTimeStampDays)),'--','Color',rgb('Black'),'LineWidth',2);
        hold off
        
        % Fit: Linear
        fitType = fittype('a*x + b'); % The equation for your fit goes here
        [fLin,gofLin] = fit(StartTimeStampDays',CrateEquivalentmV{j,i}',fitType,...
            'Weights',(1./rateEquivalentmV_E{j,i}).^2,...
            'StartPoint',[0 CrateEquivalentmV{j,i}(1)],...
            'Robust','Bisquare');
        ci = confint(fLin,0.68); uncertaintyLin = (ci(2,1)-ci(1,1))/2;
        hold on
        hl = plot(StartTimeStampDays,fLin.a*StartTimeStampDays+fLin.b,'-','Color',rgb('DarkBlue'),'LineWidth',2);
        hold off
        
        % SAVE SLOPES & FITS
        SlopeRW123PSR1234_mV{j,i} = fLin.a;
        SlopeErrorRW123PSR1234_mV{j,i} = uncertaintyLin;
        FitsRW123PSR1234{j,i} = fLin.a*StartTimeStampDays+fLin.b;
        
        % Fit: Poly3
        fitType = fittype('a*x^3 + b*x^2 + c*x + d'); % The equation for your fit goes here
        [f,gofQuad] = fit(StartTimeStampDays',CrateEquivalentmV{j,i}',fitType,...
            'Weights',(1./rateEquivalentmV_E{j,i}).^2,...
            'StartPoint',[0 0 0 CrateEquivalentmV{j,i}(1)],...
            'Robust','Bisquare');
        ci = confint(f,0.68); uncertaintyQuad = (ci(2,2)-ci(1,2))/2;
        hold on
        hq = plot(StartTimeStampDays,...
            f.a*StartTimeStampDays.^3 + f.b*StartTimeStampDays.^2 + f.c*StartTimeStampDays + f.d,...
            '-','Color',rgb('DarkGreen'),'LineWidth',2);
        hold off
        
        ylabel('\Delta U (meV)');
        xlabel('Days since last RW setting');
        leg=legend([hi hc hl hq],...
            sprintf('Pseudo-Ring %0.f - std: %.2g meV ',i,std(CrateEquivalentmV{j,i})),...
            sprintf('constant - p-val=%.2g',chi2pvalue(gofCons.sse,gofCons.dfe)),...
            sprintf('linear: %.1f+-%.1f mV/day - p-val=%.2g',fLin.a,uncertaintyLin,chi2pvalue(gofLin.sse,gofLin.dfe)),...
            sprintf('poly3 - p-val=%.2g',chi2pvalue(gofQuad.sse,gofQuad.dfe)),...
            'location','eastoutside','FontSize',12);
        leg.Color = 'none'; legend boxoff;
        PrettyFigureFormat
    end
    
    %% Rate Evolution
    myMainTitle = sprintf('KATRIN - KNM2 RW%s - FPD Rate Evolution 300eV below Endpoint',num2str(j));
    maintitle   = myMainTitle;
    savefile2    = sprintf('plots/KNM2_RM300_EffectivePotentialFluctuation_RW%.0f_PseudoRings_2.png',j);
    fig2      = figure('Name',sprintf('KATRIN - %s Scanwise Background',RunList),...
        'NumberTitle','off','rend','painters','pos',[10 10 1400 1000]);
    a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
    a.FontSize=24;a.FontWeight='bold';
    
    for i=1:A.nRings
        
        subplot(A.nRings,1,i)
        hnc=plot(R.SingleRunData.StartTimeStamp,rateRaw(i,:),'s--','Color',rgb('Green'),'LineWidth',1,'MarkerSize',8,'markerfacecolor',rgb('DarkGreen'));
        hold on
        hc=plot(R.SingleRunData.StartTimeStamp,rate(i,:),'s--','Color',rgb('Red'),'LineWidth',1,'MarkerSize',12,'markerfacecolor',rgb('IndianRed'));
        hold off
        ylabel('cps');
        xlabel('Scan Start Time');
        leg=legend([hnc hc],...
            sprintf('before correction'),...
            sprintf('after correction'),...
            'location','northeast');
        leg.Color = 'none'; legend boxoff;
        PrettyFigureFormat
    end
    
    %export_fig(fig1,savefile1);
    %export_fig(fig2,savefile2);
end

%% Overall Diagram
% OverallStartTimeStamp
% CrateEquivalentmV
% rateEquivalentmV_E

  myMainTitle = sprintf('KATRIN KNM2 - FPD Rate 300eV below Endpoint - Average WGTS Potential Drift');
    maintitle   = myMainTitle;
    savefile3    = sprintf('plots/KNM2_RM300_EffectivePotentialFluctuation_RW%.0f_PseudoRings_3.png',123);
    fig3      = figure('Name',sprintf('KATRIN - %s Scanwise Background',RunList),...
        'NumberTitle','off','rend','painters','pos',[10 10 1400 1000]);
    a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
    a.FontSize=24;a.FontWeight='bold';
    
    
    % ABSOLUTE SHIFTS CORRECTION
    
for i=1:A.nRings
    TimeRW123        = [OverallStartTimeStamp{1} OverallStartTimeStamp{2} OverallStartTimeStamp{3}];
    RatemVRW123      = [CrateEquivalentmV{1,i} CrateEquivalentmV{2,i} CrateEquivalentmV{3,i}];
    RateErrormVRW123 = [rateEquivalentmV_E{1,i} rateEquivalentmV_E{2,i} rateEquivalentmV_E{3,i}];
    FitsRW123        = [FitsRW123PSR1234{1,i} FitsRW123PSR1234{2,i} FitsRW123PSR1234{3,i}];
        
    subplot(A.nRings,1,i)
    H1=plot(TimeRW123,RatemVRW123,'s--','Color',rgb('Green'),'LineWidth',1,'MarkerSize',8,'markerfacecolor',rgb('IndianRed'));
    
    % Fit: Linear
    fitType = fittype('a*x + b'); % The equation for your fit goes here
    [fLin,gofLin] = fit((days(TimeRW123-TimeRW123(1)))',RatemVRW123',fitType,...
        'Weights',(1./RateErrormVRW123').^2,...
        'StartPoint',[0 RatemVRW123(1)],...
        'Robust','Bisquare');
    ci = confint(fLin,0.68); uncertaintyLin = (ci(2,1)-ci(1,1))/2;
    hold on
    H2 = plot(OverallStartTimeStamp{1} ,FitsRW123PSR1234{1,i},'-','Color',rgb('DarkBlue'),'LineWidth',2);
    H3 = plot(OverallStartTimeStamp{2} ,FitsRW123PSR1234{2,i},'-','Color',rgb('DarkRed'),'LineWidth',2);
    H4 = plot(OverallStartTimeStamp{3} ,FitsRW123PSR1234{3,i},'-','Color',rgb('DarkOrange'),'LineWidth',2);
    hold off
    
    ylabel('\Delta U (meV)');
        xlabel('Days');
        leg=legend([H1 H2 H3 H4],...
            sprintf('Pseudo-Ring %0.f',i),...
            sprintf('linear RW1: %.1f+-%.1f mV/day',SlopeRW123PSR1234_mV{1,i},SlopeErrorRW123PSR1234_mV{1,i}),...
            sprintf('linear RW2: %.1f+-%.1f mV/day',SlopeRW123PSR1234_mV{2,i},SlopeErrorRW123PSR1234_mV{2,i}),...
            sprintf('linear RW3: %.1f+-%.1f mV/day',SlopeRW123PSR1234_mV{3,i},SlopeErrorRW123PSR1234_mV{3,i}),...
            'location','eastoutside','FontSize',12);
            
            
        leg.Color = 'none'; legend boxoff;
        PrettyFigureFormat
end

save('SamakKNM2_DriftInRW123PSR1234_mVperDay.mat','SlopeRW123PSR1234_mV','SlopeErrorRW123PSR1234_mV');

%% Correction Difference between methods
% OverallStartTimeStamp
% CrateEquivalentmV
% rateEquivalentmV_E

  myMainTitle = sprintf('KATRIN KNM2 - Correction Difference');
    maintitle   = myMainTitle;
    savefile3    = sprintf('plots/KNM2_RM300_CorrectionDiff_RW%.0f_PseudoRings_3.png',123);
    fig4      = figure('Name',sprintf('KATRIN - %s Correction Diff',RunList),...
        'NumberTitle','off','rend','painters','pos',[10 10 1400 1000]);
    a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
    a.FontSize=24;a.FontWeight='bold';
    
    
    % ABSOLUTE SHIFTS CORRECTION
    
for i=1:A.nRings
    TimeRW123        = [OverallStartTimeStamp{1} OverallStartTimeStamp{2} OverallStartTimeStamp{3}];
    Rate1RW123       = [Crate1{1,i} Crate1{2,i} Crate1{3,i}];
    Rate2RW123       = [Crate2{1,i} Crate2{2,i} Crate2{3,i}];
    percentages      = std(Rate2RW123-Rate1RW123)./mean(rateRaw(:,i));
        
    subplot(A.nRings,1,i)
    H1=plot(TimeRW123,(Rate2RW123-Rate1RW123),'s--','Color',rgb('Green'),'LineWidth',1,'MarkerSize',8,'markerfacecolor',rgb('IndianRed'));
    
    ylabel('(\Delta r)/r (%)');
    ylim([-inf inf]);
        xlabel('Days');
        leg=legend([H1],...
            sprintf('Pseudo-Ring %0.f',i),'location','southeast','FontSize',12);
            
        leg.Color = 'none'; legend boxoff;
        PrettyFigureFormat
end

%% Overall Diagram
% OverallStartTimeStamp
% CrateEquivalentmV
% rateEquivalentmV_E

  myMainTitle = sprintf('KATRIN KNM2 - FPD Rate E_0-300eV -  WGTS Potential Drift Versus Cumulative Throughput');
    maintitle   = myMainTitle;
    savefile3    = sprintf('plots/KNM2_RM300_EffectivePotentialFluctuation_RW%.0f_PseudoRings_3.png',123);
    fig5      = figure('Name',sprintf('KATRIN - %s Scanwise Background',RunList),...
        'NumberTitle','off','rend','painters','pos',[10 10 1400 1000]);
    a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
    a.FontSize=24;a.FontWeight='bold';
    
    
for i=1:A.nRings
    TimeRW123        = [OverallStartTimeStamp{1} OverallStartTimeStamp{2} OverallStartTimeStamp{3}];
    RatemVRW123      = [CrateEquivalentmV{1,i} CrateEquivalentmV{2,i} CrateEquivalentmV{3,i}];
    RateErrormVRW123 = [rateEquivalentmV_E{1,i} rateEquivalentmV_E{2,i} rateEquivalentmV_E{3,i}];
        
    CumTActivity     = cumsum([0 diff(days(TimeRW123-TimeRW123(1)))] .* [Activity{1} Activity{2} Activity{3}]);
    CumTActivity     = CumTActivity./5e17;
    
    subplot(A.nRings,1,i)
    H1=plot(CumTActivity,RatemVRW123,'s--','Color',rgb('Green'),'LineWidth',1,'MarkerSize',8,'markerfacecolor',rgb('IndianRed'));
    
    % Fit: Linear
    fitType = fittype('a*x + b'); % The equation for your fit goes here
    [fLin,gofLin] = fit(CumTActivity',RatemVRW123',fitType,...
        'Weights',(1./RateErrormVRW123').^2,...
        'StartPoint',[0 RatemVRW123(1)],...
        'Robust','Bisquare');
    ci = confint(fLin,0.68); uncertaintyLin = (ci(2,1)-ci(1,1))/2;
    hold on
    H2 = plot(CumTActivity,fLin.a*(CumTActivity)+fLin.b,'-','Color',rgb('DarkBlue'),'LineWidth',2);
    hold off
    
    xlabel('Cumulative Throughput Unit (day at nominal column density and T-purity)');
    ylabel('\Delta U \rho (meV)');
        leg=legend([H1 H2],...
            sprintf('Pseudo-Ring %0.f',i),...
            sprintf('linear: %.1f+-%.1f mV/CTU',fLin.a,uncertaintyLin),...
            'location','southeast','FontSize',12);
        leg.Color = 'none'; legend boxoff;
        PrettyFigureFormat
end

%% Mean Values Per Period Per Ring
% CrateEquivalentmVAverage
% i=psr
% j=period
NbxRunsPeriod = [121 95 92];
for j=1:3
    for i=1:4
        fprintf('\n Period %0.f PSR%0.f = %.1f mV\n',j,i,mean(CrateEquivalentmV{j,i}));
        CrateEquivalentmVAverage(i,j)  = mean(CrateEquivalentmV{j,i});
        CrateEquivalentmVAverageE(i,j) = mean(rateEquivalentmV_E{j,i})./sqrt(NbxRunsPeriod(j));
        SlopeEquivalent_mV(i,j)       = SlopeRW123PSR1234_mV{j,i};
        SlopeErrorEquivalent_mV(i,j)  = SlopeErrorRW123PSR1234_mV{j,i};
    end
end

%% Plot Rate Per Period
% Create a ribbon point using the ribbon function
figure
hh = ribboncoloredZ(CrateEquivalentmVAverage');
hcb=colorbar; colormap(cool);ylabel(hcb, 'mV');
hh(1).LineWidth = 3;hh(2).LineWidth = 3;hh(3).LineWidth = 3;hh(4).LineWidth = 3;
hh(1).MeshStyle = 'both';hh(2).MeshStyle = 'both';hh(3).MeshStyle = 'both';hh(4).MeshStyle = 'both';
hh(1).Marker = '.';hh(2).Marker = '.';hh(3).Marker = '.';hh(4).Marker = '.';
hh(1).MarkerSize = 20;hh(2).MarkerSize = 20;hh(3).MarkerSize = 20;hh(4).MarkerSize = 20;
% Add title and axis labels
ht=title(sprintf('Mean Potential Per Period Per Pseudo-ring \n ROI %s - MOS %s - KNM1 Calibration %s',ROI,HVdriftCorFlag,KNM1CorFlag));
ylabel('RW Period'); ylim([0.5 3.5]) ; yticks([1 2 3]); %ylabel({'RW1';'RW2';'RW3'});
xlabel('Pseudo-Ring'); xlim([0.5 4.5]); xticks([1 2 3 4]);%xlabel({'PSR1','PSR2','PSR3','PSR4'});
zlabel('mV-equivalent'); %zlim([-102 1200]);
PrettyFigureFormat
ht.FontSize=16;
%view(90,0); % period on x-axis
%view(0,90); % period on x-axis
%view(0,-90); % period on x-axis

%% Plot Slope Rate Per Period
% Create a ribbon point using the ribbon function
figure
hh = ribboncoloredZ(SlopeEquivalent_mV');
hcb=colorbar; colormap(cool);ylabel(hcb, 'mV/day');
hh(1).LineWidth = 3;hh(2).LineWidth = 3;hh(3).LineWidth = 3;hh(4).LineWidth = 3;
hh(1).MeshStyle = 'both';hh(2).MeshStyle = 'both';hh(3).MeshStyle = 'both';hh(4).MeshStyle = 'both';
hh(1).Marker = '.';hh(2).Marker = '.';hh(3).Marker = '.';hh(4).Marker = '.';
hh(1).MarkerSize = 20;hh(2).MarkerSize = 20;hh(3).MarkerSize = 20;hh(4).MarkerSize = 20;
% Add title and axis labels
ht=title(sprintf('Slope Per Period Per Pseudo-ring \n ROI %s - MOS %s - KNM1 Calibration %s',ROI,HVdriftCorFlag,KNM1CorFlag));
ylabel('RW Period'); ylim([0.5 3.5]) ; yticks([1 2 3]); %ylabel({'RW1';'RW2';'RW3'});
xlabel('Pseudo-Ring'); xlim([0.5 4.5]); xticks([1 2 3 4]);%xlabel({'PSR1','PSR2','PSR3','PSR4'});
zlabel('mV-equivalent / day'); %zlim([0 10]);
PrettyFigureFormat
ht.FontSize=16;