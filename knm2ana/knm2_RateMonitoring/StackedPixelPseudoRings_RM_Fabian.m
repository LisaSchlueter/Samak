KNM1CorFlag    = 'ON';
HVdriftCorFlag = 'OFF';
ROI            = 'Default';     %Default, 14keV
Corr           = 'Fabian';      %Fabian, Thierry
CorrMode       = 'MC';         %Fit, MC
Mode           = 'Periodwise';    %Ringwise, Periodwise

%% KNM1 with Calibration - Divide Rates by
switch KNM1CorFlag
    case 'OFF'
        KNM1correction  = [ 1.0000    1         1         1     ];
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

if strcmp(Mode,'Ringwise')
    %
    % Reference Rate: KNM2_RW2
    %
        DataType  = 'Real';
        RunAnaArg = {'RunList','KNM2_RW2','DataType',DataType,...
            'FSDFlag','BlindingKNM2','ELossFlag','KatrinT2',...
            'AnaFlag','StackPixel','RingMerge','Full'};
        R        = MultiRunAnalysis(RunAnaArg{:});
        A         = RingAnalysis('RunAnaObj',R,'RingList',1:4);
        R         = A.MultiObj(1);
        R.ROIFlag=ROI; R.SetROI;


        %% Time in days
        %StartTimeStampDays           = days(R.SingleRunData.StartTimeStamp-R.SingleRunData.StartTimeStamp(1));

        %% HV Drift Correction - MOS
        FirstDayPeriod1 = datetime('02-Oct-2019 14:52:19');
        TimeLineDaysFirstDayPeriod1 = days(R.SingleRunData.StartTimeStamp-FirstDayPeriod1);
        HVdriftPerPixel = HVdriftMvDay*(TimeLineDaysFirstDayPeriod1) * 6.303e-3;

        %% Stacked Pixel Data for each patch
        count  = zeros(A.nRings,numel(A.RunAnaObj.RunList));
        rateRef   = zeros(A.nRings,numel(A.RunAnaObj.RunList));
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
            R.RMCorrection('saveplot','OFF','pixlist',sprintf('ring%i',i),'QAplots','OFF','Mode',CorrMode);

            % Include HV Drift if any
            count(i,:)  = R.SingleRunData.TBDIS_RM + HVdriftPerPixel.*numel(R.PixList).*R.SingleRunData.qUfrac_RM.*R.SingleRunData.TimeSec;

            % Correct for KNM1 Radial Effect if any
            count(i,:)  = count(i,:) ./ KNM1correction(i);

            sstime(i,:) = R.SingleRunData.qUfrac_RM.*R.SingleRunData.TimeSec;
            rateRef(i,:)   = count(i,:)./sstime(i,:);
            rateE(i,:)  = sqrt(count(i,:))./sstime(i,:);
            RefPeriodRW2CPS(i) = mean(rateRef(i,:));
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

%         if j==1
%             RunList = 'KNM2_RW1a';
%         elseif j==2
%             RunList = 'KNM2_RW1b';
%         else
%             RunList   = ['KNM2_RW' num2str(j-1)];
%         end
        RunList   = ['KNM2_RW' num2str(j)];

        % Read Data
        DataType  = 'Real';
        RunAnaArg = {'RunList',RunList,'DataType',DataType,...
            'FSDFlag','BlindingKNM2','ELossFlag','KatrinT2',...
            'NonPoissonScaleFactor',1,'AnaFlag','StackPixel','RingMerge','Full'};
        R        = MultiRunAnalysis(RunAnaArg{:});
        A         = RingAnalysis('RunAnaObj',R,'RingList',1:4);
        R         = A.MultiObj(1);
        R.ROIFlag=ROI; R.SetROI;

        % Time in days
        StartTimeStampDays{j}        = days(R.SingleRunData.StartTimeStamp-R.SingleRunData.StartTimeStamp(1));
        OverallStartTimeStamp{j}     = (R.SingleRunData.StartTimeStamp);

        % HV Drift Correction
        FirstDayPeriod1 = datetime('02-Oct-2019 14:52:19');
        TimeLineDaysFirstDayPeriod1 = days(R.SingleRunData.StartTimeStamp-FirstDayPeriod1);
        HVdriftPerPixel  = HVdriftMvDay*(TimeLineDaysFirstDayPeriod1) * 6.303e-3;

        % Activity
        Activity{j} =  (R.SingleRunData.WGTS_MolFrac_TT+0.5*R.SingleRunData.WGTS_MolFrac_HT+0.5*R.SingleRunData.WGTS_MolFrac_DT)...
            .*(R.SingleRunData.WGTS_CD_MolPerCm2);

        % Stacked Pixel Data for each patch
        count   = zeros(A.nRings,numel(A.RunAnaObj.RunList));
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
            rateRaw{j,i} = count(i,:)./sstime(i,:);

            % Ref = Period2 - Average qU correction  
            HVcorrCPSperPixel{j,i} = (qUmeanRW{j,i} - qUmeanRW2(i)) * 6.3032 * numel(R.PixList);

            % Correct for activity and qU 1
            cf(i,:)     = R.RMRateErosCorrectionqUActivity;
            Crate1{j,i} = (rateRaw{j,i}.*cf(i,:)./mean(Activity{j}).*ActivityRW2 + HVcorrCPSperPixel{j,i});% - (RefPeriodRW2CPS(1)./numel(A.MultiObj(1).PixList).*numel(R.PixList));
            

            % Correct for activity and qU
            R.RMCorrection('saveplot','ON','pixlist',sprintf('ring%i',i),'QAplots','OFF','Mode',CorrMode);                                                 % KNM1 correction included here!
            count(i,:)  = R.SingleRunData.TBDIS_RM + HVdriftPerPixel.*numel(R.PixList).*R.SingleRunData.qUfrac_RM.*R.SingleRunData.TimeSec;

            rate{j,i}   = count(i,:)./sstime(i,:);
            rateE(i,:)  = sqrt(count(i,:))./sstime(i,:);



            % Ref = Period2 PSR1
            % j = RW Period
            % i = PSR
            Crate2{j,i} =  (rate{j,i}./mean(Activity{j}).*ActivityRW2 + HVcorrCPSperPixel{j,i});% - (RefPeriodRW2CPS(1)./numel(A.MultiObj(1).PixList).*numel(R.PixList));
            
            if strcmp(Corr,'Fabian')
                CrateEquivalentmV{j,i}    =  -(Crate2{j,i}) ./737.8 * 1e3 * 117 / numel(R.PixList);
            elseif strcmp(Corr,'Thierry')
                CrateEquivalentmV{j,i}    =  -(Crate1{j,i}) ./737.8 * 1e3 * 117 / numel(R.PixList);
            end
            
            rateEquivalentmV_E{j,i}   =  (rateE(i,:) ./737.8 *1e3 * 117 / numel(R.PixList));
            
            
            % Fit: Constant
            fitType = fittype('a*x+b'); % The equation for your fit goes here
            [f,gofCons] = fit(StartTimeStampDays{j}',CrateEquivalentmV{j,i}',fitType,...
                'Weights',(1./rateEquivalentmV_E{j,i}).^2,...
                'StartPoint',[0 CrateEquivalentmV{j,i}(1)],...
                'Robust','Bisquare',...
                'Lower',[0 -1e9],...
                'Upper',[0 1e9]);
            ci = confint(f,0.68); uncertaintyCons = (ci(2,1)-ci(1,1))/2;

            % Fit: Linear
            fitType = fittype('a*x + b'); % The equation for your fit goes here
            [fLin,gofLin] = fit(StartTimeStampDays{j}',CrateEquivalentmV{j,i}',fitType,...
                'Weights',(1./rateEquivalentmV_E{j,i}).^2,...
                'StartPoint',[0 CrateEquivalentmV{j,i}(1)],...
                'Robust','Bisquare');
            ci = confint(fLin,0.68); uncertaintyLin = (ci(2,1)-ci(1,1))/2;

            % Fit: Poly3
            fitType = fittype('a*x^3 + b*x^2 + c*x + d'); % The equation for your fit goes here
            [f,gofQuad] = fit(StartTimeStampDays{j}',CrateEquivalentmV{j,i}',fitType,...
                'Weights',(1./rateEquivalentmV_E{j,i}).^2,...
                'StartPoint',[0 0 0 CrateEquivalentmV{j,i}(1)],...
                'Robust','Bisquare');
            ci = confint(f,0.68); uncertaintyQuad = (ci(2,2)-ci(1,2))/2;
            
            % SAVE SLOPES & FITS
            SlopeRW123PSR1234_mV{j,i} = fLin.a;
            pValCons{j,i} = chi2pvalue(gofCons.sse,gofCons.dfe);
            SlopeErrorRW123PSR1234_mV{j,i} = uncertaintyLin;
            LinRW123PSR1234{j,i} = fLin.a*StartTimeStampDays{j}+fLin.b;
            pValLin{j,i} = chi2pvalue(gofLin.sse,gofLin.dfe);
            Poly3RW123PSR1234{j,i} = f.a*StartTimeStampDays{j}.^3 + f.b*StartTimeStampDays{j}.^2 + f.c*StartTimeStampDays{j} + f.d;
            pValPoly3{j,i} = chi2pvalue(gofQuad.sse,gofQuad.dfe);

            
        end
    end
    
    %% Mean Values Per Period Per Ring
    % CrateEquivalentmVAverage
    % i=psr
    % j=period
    NbxRunsPeriod = [121 95 92];
    for j=1:3
        for i=1:4
            CrateEquivalentmVAverage(i,j)  = mean(CrateEquivalentmV{j,i});
            CrateEquivalentmVAverageE(i,j) = mean(rateEquivalentmV_E{j,i})./sqrt(NbxRunsPeriod(j));
            SlopeEquivalent_mV(i,j)       = SlopeRW123PSR1234_mV{j,i};
            SlopeErrorEquivalent_mV(i,j)  = SlopeErrorRW123PSR1234_mV{j,i};
        end
    end
    
    CrateEquivalentmVAverage = CrateEquivalentmVAverage-mean(mean(CrateEquivalentmVAverage));
    for j=1:3
        for i=1:4
            fprintf('\n Period %0.f PSR%0.f = (%.1f +- %.1f) mV\n',j,i,CrateEquivalentmVAverage(i,j),CrateEquivalentmVAverageE(i,j));
        end
    end
    
    save('SamakKNM2_CratemVequivalentInRW123PSR1234.mat','rateRaw','rate','CrateEquivalentmV','rateEquivalentmV_E','OverallStartTimeStamp','StartTimeStampDays');
    save('SamakKNM2_CratemVequivalentSlopesAndFitsInRW123PSR1234.mat','SlopeRW123PSR1234_mV','pValCons','pValLin','pValPoly3','SlopeErrorRW123PSR1234_mV','LinRW123PSR1234','Poly3RW123PSR1234');
    save('SamakKNM2_ShiftDriftInRW123PSR1234_mVperDay.mat','CrateEquivalentmVAverage','CrateEquivalentmVAverageE','SlopeEquivalent_mV','SlopeErrorEquivalent_mV');
    
end

if strcmp(Mode,'Periodwise')
    %
    % Reference Rate: KNM2_RW2
    %
        DataType  = 'Real';
        RunAnaArg = {'RunList','KNM2_RW2','DataType',DataType,...
            'FSDFlag','BlindingKNM2','ELossFlag','KatrinT2',...
            'AnaFlag','StackPixel','RingMerge','Full'};
        Ref        = MultiRunAnalysis(RunAnaArg{:});
        Ref.ROIFlag=ROI; Ref.SetROI;


        %% Time in days
        StartTimeStampDays           = days(Ref.SingleRunData.StartTimeStamp-Ref.SingleRunData.StartTimeStamp(1));

        %% HV Drift Correction - MOS
        FirstDayPeriod1 = datetime('02-Oct-2019 14:52:19');
        TimeLineDaysFirstDayPeriod1 = days(Ref.SingleRunData.StartTimeStamp-FirstDayPeriod1);
        HVdriftPerPixel = HVdriftMvDay*(TimeLineDaysFirstDayPeriod1) * 6.303e-3;


        % Slow Control Data: HV setting
        qUCorrRW2       = (Ref.SingleRunData.qU_RM);
        qUmeanRW2       = mean(qUCorrRW2);

        % Correct for activity and qU
        Ref.RMCorrection('saveplot','OFF','QAplots','OFF','Mode',CorrMode);

        % Include HV Drift if any
        count  = Ref.SingleRunData.TBDIS_RM + HVdriftPerPixel.*numel(Ref.PixList).*Ref.SingleRunData.qUfrac_RM.*Ref.SingleRunData.TimeSec;

        % Correct for KNM1 Radial Effect if any
        count  = count ./ (28/118.*KNM1correction(1)+36/117.*KNM1correction(2)+34/117.*KNM1correction(3)+19/117.*KNM1correction(4));

        sstime = Ref.SingleRunData.qUfrac_RM.*Ref.SingleRunData.TimeSec;
        rate   = count./sstime;
        rateE  = sqrt(count)./sstime;
        RefPeriodRW2CPS = mean(rate);
        ActivityRW2  =  mean(Ref.SingleRunData.WGTS_MolFrac_TT'+0.5*Ref.SingleRunData.WGTS_MolFrac_HT'+0.5*Ref.SingleRunData.WGTS_MolFrac_DT')...
            .*mean(Ref.SingleRunData.WGTS_CD_MolPerCm2);

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
        R        = MultiRunAnalysis(RunAnaArg{:});
        R.ROIFlag=ROI; R.SetROI;

        % Time in days
        StartTimeStampDays           = days(R.SingleRunData.StartTimeStamp-R.SingleRunData.StartTimeStamp(1));
        OverallStartTimeStamp{j}     = (R.SingleRunData.StartTimeStamp);

        % HV Drift Correction
        FirstDayPeriod1 = datetime('02-Oct-2019 14:52:19');
        TimeLineDaysFirstDayPeriod1 = days(R.SingleRunData.StartTimeStamp-FirstDayPeriod1);
        HVdriftPerPixel  = HVdriftMvDay*(TimeLineDaysFirstDayPeriod1) * 6.303e-3;

        % Activity
        Activity{j} =  (R.SingleRunData.WGTS_MolFrac_TT+0.5*R.SingleRunData.WGTS_MolFrac_HT+0.5*R.SingleRunData.WGTS_MolFrac_DT)...
            .*(R.SingleRunData.WGTS_CD_MolPerCm2);


        % HV setting
        qUCorrRW        = (R.SingleRunData.qU_RM);
        qUmeanRW{j}        = mean(R.SingleRunData.qU_RM);

        sstime = R.SingleRunData.qUfrac_RM.*R.SingleRunData.TimeSec;

        % Include HV Drift
        count  = R.SingleRunData.TBDIS_RM + HVdriftPerPixel.*numel(R.PixList).*R.SingleRunData.qUfrac_RM.*R.SingleRunData.TimeSec;

        % Correct for KNM1 Radial Effect
        count  = count./(28/118.*KNM1correction(1)+36/117.*KNM1correction(2)+34/117.*KNM1correction(3)+19/117.*KNM1correction(4));

        % Rate uncorrected for activity and qU
        rateRaw = count./sstime;

        % Ref = Period2 - Average qU correction  
        HVcorrCPSperPixel{j} = (qUmeanRW{j} - qUmeanRW2) * 6.3032 * numel(R.PixList);

        % Correct for activity and qU 1
        cf     = R.RMRateErosCorrectionqUActivity;
        Crate1{j} = (rateRaw.*cf./mean(Activity{j}).*ActivityRW2 + HVcorrCPSperPixel{j});
        % Correct for activity and qU
        R.RMCorrection('saveplot','ON','QAplots','OFF','Mode',CorrMode);
        count  = R.SingleRunData.TBDIS_RM + HVdriftPerPixel.*numel(R.PixList).*R.SingleRunData.qUfrac_RM.*R.SingleRunData.TimeSec;
        count  = count./(28/118.*KNM1correction(1)+36/117.*KNM1correction(2)+34/117.*KNM1correction(3)+19/117.*KNM1correction(4));

        rate   = count./sstime;
        rateE  = sqrt(count)./sstime;



        % Ref = Period2 PSR1
        % j = RW Period
        Crate2{j}                =  (rate./mean(Activity{j}).*ActivityRW2 + HVcorrCPSperPixel{j});

        if strcmp(Corr,'Fabian')
            CrateEquivalentmV{j}    =  -(rate./mean(Activity{j}).*ActivityRW2 + HVcorrCPSperPixel{j})./737.8 * 1e3;
        elseif strcmp(Corr,'Thierry')
            CrateEquivalentmV{j}    =  -(rateRaw.*cf./mean(Activity{j}).*ActivityRW2 + HVcorrCPSperPixel{j})./737.8 * 1e3;
        end

        rateEquivalentmV_E{j}   =  (rateE ./737.8 *1e3);


        %% Rate Evolution --> mV equivalent
        myMainTitle = sprintf('KATRIN - KNM2 RW%s - FPD Rate Evolution 300eV below Endpoint',num2str(j));
        maintitle   = myMainTitle;
        savefile1    = sprintf('plots/KNM2_RM300_EffectivePotentialFluctuation_RW%.0f_PseudoRings_1.png',j);
        fig1      = figure('Name',sprintf('KATRIN - %s Scanwise Background',RunList),...
            'NumberTitle','off','rend','painters','pos',[10 10 1400 1000]);
        a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
        a.FontSize=24;a.FontWeight='bold';


        hi=plot(StartTimeStampDays,CrateEquivalentmV{j},'s--','Color',rgb('IndianRed'),'LineWidth',1,'MarkerSize',12,'markerfacecolor',rgb('IndianRed'));

        % Fit: Constant
        fitType = fittype('a*x+b'); % The equation for your fit goes here
        [f,gofCons] = fit(StartTimeStampDays',CrateEquivalentmV{j}',fitType,...
            'Weights',(1./rateEquivalentmV_E{j}).^2,...
            'StartPoint',[0 CrateEquivalentmV{1}(1)],...
            'Robust','Bisquare',...
            'Lower',[0 -1e9],...
            'Upper',[0 1e9]);
        ci = confint(f,0.68); uncertaintyCons = (ci(2,1)-ci(1,1))/2;
        hold on
        hc = plot(StartTimeStampDays,f.a.*ones(1,numel(StartTimeStampDays)),'--','Color',rgb('Black'),'LineWidth',2);
        hold off

        % Fit: Linear
        fitType = fittype('a*x + b'); % The equation for your fit goes here
        [fLin,gofLin] = fit(StartTimeStampDays',CrateEquivalentmV{j}',fitType,...
            'Weights',(1./rateEquivalentmV_E{j}).^2,...
            'StartPoint',[0 CrateEquivalentmV{j}(1)],...
            'Robust','Bisquare');
        ci = confint(fLin,0.68); uncertaintyLin = (ci(2,1)-ci(1,1))/2;
        hold on
        hl = plot(StartTimeStampDays,fLin.a*StartTimeStampDays+fLin.b,'-','Color',rgb('DarkBlue'),'LineWidth',2);
        hold off

        % SAVE SLOPES & FITS
        SlopeRW123PSR1234_mV{j} = fLin.a;
        SlopeErrorRW123PSR1234_mV{j} = uncertaintyLin;
        FitsRW123PSR1234{j} = fLin.a*StartTimeStampDays+fLin.b;

        % Fit: Poly3
        fitType = fittype('a*x^3 + b*x^2 + c*x + d'); % The equation for your fit goes here
        [f,gofQuad] = fit(StartTimeStampDays',CrateEquivalentmV{j}',fitType,...
            'Weights',(1./rateEquivalentmV_E{j}).^2,...
            'StartPoint',[0 0 0 CrateEquivalentmV{j}(1)],...
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
            sprintf('std: %.2g meV ',std(CrateEquivalentmV{j})),...
            sprintf('constant - p-val=%.2g',chi2pvalue(gofCons.sse,gofCons.dfe)),...
            sprintf('linear: %.1f+-%.1f mV/day - p-val=%.2g',fLin.a,uncertaintyLin,chi2pvalue(gofLin.sse,gofLin.dfe)),...
            sprintf('poly3 - p-val=%.2g',chi2pvalue(gofQuad.sse,gofQuad.dfe)),...
            'location','eastoutside','FontSize',12);
        leg.Color = 'none'; legend boxoff;
        PrettyFigureFormat;
        pause(1);

        %% Rate Evolution
        myMainTitle = sprintf('KATRIN - KNM2 RW%s - FPD Rate Evolution 300eV below Endpoint',num2str(j));
        maintitle   = myMainTitle;
        savefile2    = sprintf('plots/KNM2_RM300_EffectivePotentialFluctuation_RW%.0f_PseudoRings_2.png',j);
        fig2      = figure('Name',sprintf('KATRIN - %s Scanwise Background',RunList),...
            'NumberTitle','off','rend','painters','pos',[10 10 1400 1000]);
        a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
        a.FontSize=24;a.FontWeight='bold';


        hnc=plot(R.SingleRunData.StartTimeStamp,rateRaw,'s--','Color',rgb('Green'),'LineWidth',1,'MarkerSize',8,'markerfacecolor',rgb('DarkGreen'));
        hold on
        hc=plot(R.SingleRunData.StartTimeStamp,rate,'s--','Color',rgb('Red'),'LineWidth',1,'MarkerSize',12,'markerfacecolor',rgb('IndianRed'));
        hold off
        ylabel('cps');
        xlabel('Scan Start Time');
        leg=legend([hnc hc],...
            sprintf('before correction'),...
            sprintf('after correction'),...
            'location','northeast');
        leg.Color = 'none'; legend boxoff;
        PrettyFigureFormat;
        pause(1);

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

    TimeRW123        = [OverallStartTimeStamp{1} OverallStartTimeStamp{2} OverallStartTimeStamp{3}];
    RatemVRW123      = [CrateEquivalentmV{1} CrateEquivalentmV{2} CrateEquivalentmV{3}];
    RateErrormVRW123 = [rateEquivalentmV_E{1} rateEquivalentmV_E{2} rateEquivalentmV_E{3}];
    FitsRW123        = [FitsRW123PSR1234{1} FitsRW123PSR1234{2} FitsRW123PSR1234{3}];

    H1=plot(TimeRW123,RatemVRW123,'s--','Color',rgb('Green'),'LineWidth',1,'MarkerSize',8,'markerfacecolor',rgb('IndianRed'));

    % Fit: Linear
    fitType = fittype('a*x + b'); % The equation for your fit goes here
    [fLin,gofLin] = fit((days(TimeRW123-TimeRW123(1)))',RatemVRW123',fitType,...
        'Weights',(1./RateErrormVRW123').^2,...
        'StartPoint',[0 RatemVRW123(1)],...
        'Robust','Bisquare');
    ci = confint(fLin,0.68); uncertaintyLin = (ci(2,1)-ci(1,1))/2;
    hold on
    H2 = plot(OverallStartTimeStamp{1} ,FitsRW123PSR1234{1},'-','Color',rgb('DarkBlue'),'LineWidth',2);
    H3 = plot(OverallStartTimeStamp{2} ,FitsRW123PSR1234{2},'-','Color',rgb('DarkRed'),'LineWidth',2);
    H4 = plot(OverallStartTimeStamp{3} ,FitsRW123PSR1234{3},'-','Color',rgb('DarkOrange'),'LineWidth',2);
    hold off

    ylabel('\Delta U (meV)');
        xlabel('Days');
        leg=legend([H1 H2 H3 H4],...
            sprintf('Pseudo-Ring %0.f',i),...
            sprintf('linear RW1: %.1f+-%.1f mV/day',SlopeRW123PSR1234_mV{1},SlopeErrorRW123PSR1234_mV{1}),...
            sprintf('linear RW2: %.1f+-%.1f mV/day',SlopeRW123PSR1234_mV{2},SlopeErrorRW123PSR1234_mV{2}),...
            sprintf('linear RW3: %.1f+-%.1f mV/day',SlopeRW123PSR1234_mV{3},SlopeErrorRW123PSR1234_mV{3}),...
            'location','eastoutside','FontSize',12);


        leg.Color = 'none'; legend boxoff;
        PrettyFigureFormat;
        pause(1);

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

        TimeRW123        = [OverallStartTimeStamp{1} OverallStartTimeStamp{2} OverallStartTimeStamp{3}];
        Rate1RW123       = [Crate1{1} Crate1{2} Crate1{3}];
        Rate2RW123       = [Crate2{1} Crate2{2} Crate2{3}];
        percentages      = std(Rate2RW123-Rate1RW123)./mean(rateRaw);

        H1=plot(TimeRW123,(Rate2RW123-Rate1RW123),'s--','Color',rgb('Green'),'LineWidth',1,'MarkerSize',8,'markerfacecolor',rgb('IndianRed'));

        ylabel('(\Delta r)/r (%)');
        ylim([-inf inf]);
            xlabel('Days');
            leg=legend([H1],...
                sprintf('Pseudo-Ring %0.f',i),'location','southeast','FontSize',12);

            leg.Color = 'none'; legend boxoff;
            PrettyFigureFormat;
            pause(1);
            
            
    %% Mean Values Per Period Per Ring
    % CrateEquivalentmVAverage
    % i=psr
    % j=period
    NbxRunsPeriod = [121 95 92];
    for j=1:3
            CrateEquivalentmVAverage(j)  = mean(CrateEquivalentmV{j});
            CrateEquivalentmVAverageE(j) = mean(rateEquivalentmV_E{j})./sqrt(NbxRunsPeriod(j));
            SlopeEquivalent_mV(j)       = SlopeRW123PSR1234_mV{j};
            SlopeErrorEquivalent_mV(j)  = SlopeErrorRW123PSR1234_mV{j};
    end
    
    CrateEquivalentmVAverage = CrateEquivalentmVAverage-mean(CrateEquivalentmVAverage);
    for j=1:3
        fprintf('\n Period %0.f uniform = (%.1f +- %.1f) mV\n',j,CrateEquivalentmVAverage(j),CrateEquivalentmVAverageE(j));
    end
    
    save('SamakKNM2_ShiftDriftInRW123uniform_mVperDay.mat','CrateEquivalentmVAverage','CrateEquivalentmVAverageE','SlopeEquivalent_mV','SlopeErrorEquivalent_mV');

end