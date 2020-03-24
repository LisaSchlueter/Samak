% KNM2 - -300V Rate Monitor FPD
% Stakced-pixel Per Pseudo-Ring
% Evolution of Rate
% Rate Correction
% Conversion of rate in Potential Fluctuation
%
% Last Modified: 04/02/2020
% T. Lasserre
%

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

% Enbdpoint Fits
E0RWPSR   = zeros(4,3);
E0RWPSRE  = zeros(4,3);

%
% Reference Rate: KNM2_RW2
%
DataType  = 'Real';
RunAnaArg = {'RunList','KNM2_RW2','DataType',DataType,...
    'FSDFlag','BlindingKNM2','ELossFlag','KatrinT2',...
    'AnaFlag','StackPixel','RingMerge','Full',...
    'chi2','chi2Stat',...              % statistics only
    'NonPoissonScaleFactor',1.12,...
    'MosCorrFlag','OFF',...
    'ROIFlag','14keV',...    
    'DopplerEffectFlag','FSD'};

    MR        = MultiRunAnalysis(RunAnaArg{:});
    A         = RingAnalysis('RunAnaObj',MR,'RingList',1:4);
    R         = A.MultiObj(1);
    R.ROIFlag = ROI; R.SetROI;
        
    %% Time in days
    StartTimeStampDays           = days(R.SingleRunData.StartTimeStamp-R.SingleRunData.StartTimeStamp(1));
    
    %% HV Drift Correction - MOS
    %FirstDayPeriod1 = datetime('02-Oct-2019 14:52:19'); % old golden run list
    FirstDayPeriod1 = datetime('27-Sep-2019 13:32:58');  % new golden run list
    TimeLineDaysFirstDayPeriod1 = days(MR.SingleRunData.StartTimeStamp-FirstDayPeriod1);
    HVdriftPerPixel = HVdriftMvDay*(TimeLineDaysFirstDayPeriod1) * 6.303e-3;

    %% Stacked Pixel Data for each patch
    count  = zeros(A.nRings,numel(A.RunAnaObj.RunList));
    rate   = zeros(A.nRings,numel(A.RunAnaObj.RunList));
    rateE  = zeros(A.nRings,numel(A.RunAnaObj.RunList));
    cf     = zeros(A.nRings,numel(A.RunAnaObj.RunList));
    sstime = zeros(A.nRings,numel(A.RunAnaObj.RunList));
    RefPeriodRW2CPS = zeros(A.nRings,1);
    qUmeanRW2       = zeros(A.nRings,1);
    
    for i=1:A.nRings
        R           = A.MultiObj(i);
        R.ROIFlag   = ROI; R.SetROI;
        
        % Slow Control Data: HV setting
        qUCorrRW2       = R.SingleRunData.qU_RM(1,:);
        qUmeanRW2(i)    = mean(qUCorrRW2);
    
        % Include HV Drift if any
        count(i,:)  = R.SingleRunData.TBDIS_RM + ...
                      HVdriftPerPixel.*numel(R.PixList).*R.SingleRunData.qUfrac_RM.*R.SingleRunData.TimeSec;
        
        % Correct for KNM1 Radial Effect if any
        count(i,:)  = count(i,:) ./ KNM1correction(i);
        
        sstime(i,:)        = R.SingleRunData.qUfrac_RM.*R.SingleRunData.TimeSec;
        rate(i,:)          = count(i,:)./sstime(i,:);
        rateE(i,:)         = sqrt(count(i,:))./sstime(i,:);
        cf(i,:)            = R.RMRateErosCorrectionqUActivity; %qU/Activity Per Ring --> No Need for qU additional
        RefPeriodRW2CPS(i) = mean(rate(i,:).*cf(i,:));
    end
    
    ActivityRW2  =  mean(R.SingleRunData.WGTS_MolFrac_TT'+0.5*R.SingleRunData.WGTS_MolFrac_HT'+0.5*R.SingleRunData.WGTS_MolFrac_DT') ...
        .*mean(R.SingleRunData.WGTS_CD_MolPerCm2);

    %
    % End of Reference Rate: KNM2_RW2
    %
    
    %
    % Loop on Period RW1 RW2 RW3
    %
for j=1:3
    
    RunList   = ['KNM2_RW' num2str(j)];
    
    % Read Data
    DataType  = 'Real';
    RunAnaArg = {'RunList',RunList,'DataType',DataType,...
        'FSDFlag','BlindingKNM2','ELossFlag','KatrinT2',...
        'AnaFlag','StackPixel','RingMerge','Full',...
        'chi2','chi2Stat',...              % statistics only
        'NonPoissonScaleFactor',1.12,...
        'MosCorrFlag','OFF',... 
        'ROIFlag','14keV',...   
        'DopplerEffectFlag','FSD'}; 
    MR        = MultiRunAnalysis(RunAnaArg{:});
    A         = RingAnalysis('RunAnaObj',MR,'RingList',1:4);
    R         = A.MultiObj(1);
    R.ROIFlag = ROI; R.SetROI;
    
    % Time in days
    StartTimeStampDays           = days(R.SingleRunData.StartTimeStamp-R.SingleRunData.StartTimeStamp(1));
    OverallStartTimeStamp{j}     = R.SingleRunData.StartTimeStamp;
    
    % HV Drift Correction
    % FirstDayPeriod1 = datetime('02-Oct-2019 14:52:19');
    FirstDayPeriod1 = datetime('27-Sep-2019 13:32:58'); % new golden run list
    TimeLineDaysFirstDayPeriod1 = days(MR.SingleRunData.StartTimeStamp-FirstDayPeriod1);
    HVdriftPerPixel             = HVdriftMvDay*(TimeLineDaysFirstDayPeriod1) * 6.303e-3;
    
    % Activity
    Activity{j} =  (R.SingleRunData.WGTS_MolFrac_TT+0.5*R.SingleRunData.WGTS_MolFrac_HT+0.5*R.SingleRunData.WGTS_MolFrac_DT)...
        .*(R.SingleRunData.WGTS_CD_MolPerCm2);
    
    % Stacked Pixel Data for each patch
    count  = zeros(A.nRings,numel(A.RunAnaObj.RunList));
    rate   = zeros(A.nRings,numel(A.RunAnaObj.RunList));
    rateE  = zeros(A.nRings,numel(A.RunAnaObj.RunList));
    cf     = zeros(A.nRings,numel(A.RunAnaObj.RunList));
    sstime = zeros(A.nRings,numel(A.RunAnaObj.RunList));
    
    for i=1:A.nRings
        R           = A.MultiObj(i);
        R.ROIFlag   = ROI;
        R.SetROI;
        
        % Endpoint Fit
        R.Fit
        E0RWPSR(i,j)  = R.FitResult.par(2)+R.ModelObj.Q_i;
        E0RWPSRE(i,j) = R.FitResult.err(2);
        
        % HV setting
        qUCorrRW        = (R.SingleRunData.qU_RM(1,:));
        qUmeanRW{j,i}   = mean(R.SingleRunData.qU_RM(1,:));
    
        % Include HV Drift
        count(i,:)  = R.SingleRunData.TBDIS_RM + ...
                      HVdriftPerPixel.*numel(R.PixList).*R.SingleRunData.qUfrac_RM.*R.SingleRunData.TimeSec;
        
        % Correct for KNM1 Radial Effect
        count(i,:)  = count(i,:)./KNM1correction(i);
        sstime(i,:) = R.SingleRunData.qUfrac_RM.*R.SingleRunData.TimeSec;
        rate(i,:)   = count(i,:)./sstime(i,:);
        rateE(i,:)  = sqrt(count(i,:))./sstime(i,:);
        cf(i,:)     = R.RMRateErosCorrectionqUActivity;
        
        % Ref = Period2 - Average qU correction  
        HVcorrCPSperPixel{j,i} = (qUmeanRW{j,i} - qUmeanRW2(i)) * 6.3032 * numel(R.PixList);
        
        % Ref = Period 2
        % Crate{j,i}                =  (rate(i,:).*cf(i,:)./Activity{j}.*ActivityRW2 + HVcorrCPSperPixel{j,i} - RefPeriodRW2CPS(i));
        % CrateEquivalentmV{j,i}    =  -(rate(i,:).*cf(i,:)./mean(Activity{j}).*ActivityRW2 + HVcorrCPSperPixel{j,i} - RefPeriodRW2CPS(i)) ./737.8 * 1e3 * 117 / numel(R.PixList);
        % rateEquivalentmV_E{j,i}   =  (rateE(i,:) ./737.8 *1e3 * 117 / numel(R.PixList));

        % Ref = Period2 PSR1
        % j = RW Period
        % i = PSR
        Crate{j,i}                 =  (rate(i,:).*cf(i,:)./Activity{j}.*ActivityRW2 + HVcorrCPSperPixel{j,i} ...
                                     - RefPeriodRW2CPS(1)./numel(A.MultiObj(1).PixList).*numel(R.PixList));
        
%         CrateEquivalentmV{j,i}    =  -(rate(i,:).*cf(i,:)./mean(Activity{j}).*ActivityRW2 + HVcorrCPSperPixel{j,i} ...
%                                      - RefPeriodRW2CPS(1)./numel(A.MultiObj(1).PixList).*numel(R.PixList)) ./737.8 * 1e3 * 117 / numel(R.PixList);
%         

        CrateEquivalentmV{j,i}    =  -(rate(i,:).*cf(i,:)./mean(Activity{j}).*ActivityRW2 + ...                % Rate Renormalized to RW2 Activity
                                     HVcorrCPSperPixel{j,i} ...                                                % HV Drift Correction
                                     - RefPeriodRW2CPS(1)./numel(A.MultiObj(1).PixList).*numel(R.PixList)) ... % Reference RW2 PSR1
                                     ./ (6.3060 * numel(R.PixList)) * 1e3 ;                                    % Conversion to meV
                                 
%        rateEquivalentmV_E{j,i}   =  (rateE(i,:) ./737.8 *1e3 * 117 / numel(R.PixList));
        rateEquivalentmV_E{j,i}   =  (rateE(i,:) ./ (6.3060 * numel(R.PixList)) * 1e3);
        
    end
    
    %% Rate Evolution --> mV equivalent
    myMainTitle = sprintf('KATRIN - %s - FPD Rate Evolution 300eV below Endpoint',RunList);
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
        
        % SAVE SLOPES
        SlopeRW123PSR1234_mV{j,i} = fLin.a;
        SlopeErrorRW123PSR1234_mV{j,i} = uncertaintyLin;
        
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
    myMainTitle = sprintf('KATRIN - %s - FPD Rate Evolution 300eV below Endpoint',RunList);
    maintitle   = myMainTitle;
    savefile2    = sprintf('plots/KNM2_RM300_EffectivePotentialFluctuation_RW%.0f_PseudoRings_2.png',j);
    fig2      = figure('Name',sprintf('KATRIN - %s Scanwise Background',RunList),...
        'NumberTitle','off','rend','painters','pos',[10 10 1400 1000]);
    a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
    a.FontSize=24;a.FontWeight='bold';
    
    for i=1:A.nRings
        
        subplot(A.nRings,1,i)
        hnc=plot(R.SingleRunData.StartTimeStamp,rate(i,:),'s--','Color',rgb('Green'),'LineWidth',1,'MarkerSize',8,'markerfacecolor',rgb('DarkGreen'));
        hold on
        hc=plot(R.SingleRunData.StartTimeStamp,rate(i,:).*cf(i,:),'s--','Color',rgb('Red'),'LineWidth',1,'MarkerSize',12,'markerfacecolor',rgb('IndianRed'));
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
    
    export_fig(fig1,savefile1);
    export_fig(fig2,savefile2);
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
    H2 = plot(TimeRW123,fLin.a*(days(TimeRW123-TimeRW123(1)))+fLin.b,'-','Color',rgb('DarkBlue'),'LineWidth',2);
    hold off
    
    ylabel('\Delta U (meV)');
        xlabel('Days');
        leg=legend([H1 H2],...
            sprintf('Pseudo-Ring %0.f',i),...
            sprintf('linear: %.1f+-%.1f mV/day',fLin.a,uncertaintyLin),...
            'location','southeast','FontSize',12);
        leg.Color = 'none'; legend boxoff;
        PrettyFigureFormat
end

save('./data/SamakKNM2_DriftInRW123PSR1234_mVperDay.mat','SlopeRW123PSR1234_mV','SlopeErrorRW123PSR1234_mV');

%% Overall Diagram
% OverallStartTimeStamp
% CrateEquivalentmV
% rateEquivalentmV_E

  myMainTitle = sprintf('KATRIN KNM2 - FPD Rate E_0-300eV -  WGTS Potential Drift Verus Cumultative Throughput');
    maintitle   = myMainTitle;
    savefile3    = sprintf('plots/KNM2_RM300_EffectivePotentialFluctuation_RW%.0f_PseudoRings_3.png',123);
    fig3      = figure('Name',sprintf('KATRIN - %s Scanwise Background',RunList),...
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
NbxRunsPeriod = [171 95 92];
for j=1:3
    for i=1:4
        fprintf('\n Period %0.f PSR%0.f = %.1f mV\n',j,i,mean(CrateEquivalentmV{j,i}));
        CrateEquivalentmVAverage(i,j)  = mean(CrateEquivalentmV{j,i});
        CrateEquivalentmVAverageE(i,j) = mean(rateEquivalentmV_E{j,i})./sqrt(NbxRunsPeriod(j));
        SlopeEquivalent_mV(i,j)        = SlopeRW123PSR1234_mV{j,i};
        SlopeErrorEquivalent_mV(i,j)   = SlopeErrorRW123PSR1234_mV{j,i};
    end
end

save('./data/SamakKNM2_ShiftDriftInRW123PSR1234_mVperDay.mat','CrateEquivalentmVAverage','CrateEquivalentmVAverageE','SlopeEquivalent_mV','SlopeErrorEquivalent_mV');

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
zlabel('mV-equivalent / day'); zlim([0 10]);
PrettyFigureFormat
ht.FontSize=16;

%% Comparison with the endpoint fitss
load('./data/SamakKNM2_ShiftDriftInRW123PSR1234_mVperDay.mat');
% E0RWPSR  = 18573 + [ .681740 .742798 .572836 ;...
%     .657357 .717496 .532254 ; ...
%     .622973 .660231 .529942 ;...
%     .637982 .664792 .480520 ];
E0RWPSR_Plot =  E0RWPSR - E0RWPSR(1,2);
E0RWPSR_Plot = -E0RWPSR_Plot*1e3;

% Plot - Absolute Values - normalized to Ring 1 Period 2
myMainTitle = sprintf('KATRIN KNM2 - FPD Rate E_0-300eV Compared To Endpoint Stacked-Scan Fit');
maintitle   = myMainTitle;
savefile3   = sprintf('plots/KNM2_RM300_E0+300V_RW%.0f_4PseudoRings.png',123);
fig3        = figure('Name',sprintf('KATRIN - FPD Rate E_0-300eV - Endpoint Stacked-Scan Fit'),...
    'NumberTitle','off','rend','painters','pos',[10 10 1300 900]);
a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
a.FontSize=24;a.FontWeight='bold';
PlotStyle1 = { 'o','MarkerSize',10,'MarkerFaceColor',rgb('SkyBlue'),'LineWidth',4,'Color',rgb('SkyBlue')};
PlotStyle2 = { 's','MarkerSize',10,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',4,'Color',rgb('IndianRed')};
PlotStyle3 = { 'd','MarkerSize',10,'MarkerFaceColor',rgb('DarkBlue'),'LineWidth',4,'Color',rgb('DarkBlue')};

for i=1:3
subplot(3,1,i)
e1 = errorbar([1 2 3 4],E0RWPSR_Plot(:,i),E0RWPSRE(:,i)*1e3,PlotStyle1{:});
e1.CapSize = 0;
hold on
e2 = errorbar([1 2 3 4],CrateEquivalentmVAverage(:,i),CrateEquivalentmVAverageE(:,i)*sqrt(NbxRunsPeriod(i)),PlotStyle2{:});
e2.CapSize = 0;
hold off
xlabel(sprintf('Pseudo-Ring'));
ylabel(sprintf('RW%0.f - \\Delta U_{eq} (meV)',i));
legend([e1 e2],'Endpoint E_0+[-90;+50] eV','Monitoring Data E_0-300 eV','Location','SouthEast');
 legend boxoff;
PrettyFigureFormat;
xlim([0.5 4.5]);
%ylim([-100 200]);
xticks([1 2 3 4]);
end


% Plot - differences E0-300V - normalized to Ring 1 Period 2
myMainTitle = sprintf('KATRIN KNM2 - FPD Rate E_0-300eV Compared To Endpoint Stacked-Scan Fit');
maintitle   = myMainTitle;
savefile3   = sprintf('plots/KNM2_RM300_E0+300V_RW%.0f_4PseudoRings.png',123);
fig3        = figure('Name',sprintf('KATRIN - FPD Rate E_0-300eV - Endpoint Stacked-Scan Fit'),...
    'NumberTitle','off','rend','painters','pos',[10 10 1300 900]);
a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
a.FontSize=24;a.FontWeight='bold';
PlotStyle1 = { 'o','MarkerSize',10,'MarkerFaceColor',rgb('SkyBlue'),'LineWidth',4,'Color',rgb('SkyBlue')};
PlotStyle2 = { 's','MarkerSize',10,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',4,'Color',rgb('IndianRed')};

for i=1:3
subplot(3,1,i)
e1 = errorbar([1 2 3 4],E0RWPSR_Plot(:,i)-CrateEquivalentmVAverage(:,i),sqrt((E0RWPSRE(:,i)*1e3).^2+(CrateEquivalentmVAverageE(:,i)*sqrt(NbxRunsPeriod(i))).^2),PlotStyle3{:});
e1.CapSize = 0;
hold on
e2 = plot([1 2 3 4],E0RWPSR_Plot(:,i)*0,'LineWidth',2,'Color','Black','LineStyle','--');
hold off
xlabel(sprintf('Pseudo-Ring'));
ylabel(sprintf('RW%0.f - \\Delta U_{eq} (meV)',i));
legend([e1],'Difference Endpoint E_0+[-90;+50] eV - Monitoring Data E_0-300 eV','Location','SouthWest');
 legend boxoff;
PrettyFigureFormat;
xlim([0.5 4.5]);
%ylim([-100 200]);
xticks([1 2 3 4]);
end

    