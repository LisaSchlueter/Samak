% KNM2 - -300V Rate Monitor FPD
% Stakced-pixel Per Pseudo-Ring
% Evolution of Rate
% Rate Correction
% Conversion of rate in Potential Fluctuation
%
% Last Modified: 04/02/2020
% T. Lasserre
%

% Reference Rate
    DataType  = 'Real';
    RunAnaArg = {'RunList','KNM2_RW2','DataType',DataType,...
        'FSDFlag','BlindingKNM2','ELossFlag','KatrinT2',...
        'AnaFlag','StackPixel','RingMerge','Full'};
    MR        = MultiRunAnalysis(RunAnaArg{:});
    A         = RingAnalysis('RunAnaObj',MR,'RingList',1:4);
    R         = A.MultiObj(1);
    %% Slow Control Data
    p1 =(R.SingleRunData.WGTS_MolFrac_TT'+0.5*R.SingleRunData.WGTS_MolFrac_HT'+0.5*R.SingleRunData.WGTS_MolFrac_DT')./mean((R.SingleRunData.WGTS_MolFrac_TT'+0.5*R.SingleRunData.WGTS_MolFrac_HT'+0.5*R.SingleRunData.WGTS_MolFrac_DT')).*R.SingleRunData.WGTS_CD_MolPerCm2'./mean(R.SingleRunData.WGTS_CD_MolPerCm2');
    p2 = mean(R.SingleRunData.qU_RM,1); p2=p2-mean(p2);
    %% Time in days
    StartTimeStampDays           = days(R.SingleRunData.StartTimeStamp-R.SingleRunData.StartTimeStamp(1));
    
    %% Stacked Pixel Data for each patch
    count = zeros(A.nRings,numel(A.RunAnaObj.RunList));
    rate  = zeros(A.nRings,numel(A.RunAnaObj.RunList));
    rateE = zeros(A.nRings,numel(A.RunAnaObj.RunList));
    cf    = zeros(A.nRings,numel(A.RunAnaObj.RunList));
    RefPeriodRW2 = zeros(A.nRings,1);
    for i=1:A.nRings
        R           = A.MultiObj(i);
        count(i,:)  = R.SingleRunData.TBDIS_RM;
        sstime      = mean(R.SingleRunData.qUfrac_RM,1).*R.SingleRunData.TimeSec;
        rate(i,:)   = count(i,:)./sstime;
        rateE(i,:)  = sqrt(count(i,:))./sstime;
        cf(i,:)     = R.RMRateErosCorrectionqUActivity;
        RefPeriodRW2(i) = mean(rate(i,:).*cf(i,:));
    end
    ActivityRW2  =  mean(R.SingleRunData.WGTS_MolFrac_TT'+0.5*R.SingleRunData.WGTS_MolFrac_HT'+0.5*R.SingleRunData.WGTS_MolFrac_DT').*mean(R.SingleRunData.WGTS_CD_MolPerCm2);

  
% Loop on Period
for j=1:3
    
    RunList   = ['KNM2_RW' num2str(j)];
    
    %% Read Data
    DataType  = 'Real';
    RunAnaArg = {'RunList',RunList,'DataType',DataType,...
        'FSDFlag','BlindingKNM2','ELossFlag','KatrinT2',...
        'AnaFlag','StackPixel','RingMerge','Full'};
    MR        = MultiRunAnalysis(RunAnaArg{:});
    A         = RingAnalysis('RunAnaObj',MR,'RingList',1:4);
    R         = A.MultiObj(1);
    
    %% Slow Control Data
    p1 = (R.SingleRunData.WGTS_MolFrac_TT'+0.5*R.SingleRunData.WGTS_MolFrac_HT'+0.5*R.SingleRunData.WGTS_MolFrac_DT')./mean((R.SingleRunData.WGTS_MolFrac_TT'+0.5*R.SingleRunData.WGTS_MolFrac_HT'+0.5*R.SingleRunData.WGTS_MolFrac_DT')).*R.SingleRunData.WGTS_CD_MolPerCm2'./mean(R.SingleRunData.WGTS_CD_MolPerCm2');
    p2 = mean(R.SingleRunData.qU_RM,1); p2=p2-mean(p2);
    
    %% Time in days
    StartTimeStampDays           = days(R.SingleRunData.StartTimeStamp-R.SingleRunData.StartTimeStamp(1));
    OverallStartTimeStamp{j}     = (R.SingleRunData.StartTimeStamp);
    
    %% Activity
    Activity{j}                  =  (R.SingleRunData.WGTS_MolFrac_TT+0.5*R.SingleRunData.WGTS_MolFrac_HT+0.5*R.SingleRunData.WGTS_MolFrac_DT).*(R.SingleRunData.WGTS_CD_MolPerCm2);

    %% Stacked Pixel Data for each patch
    count = zeros(A.nRings,numel(A.RunAnaObj.RunList));
    rate  = zeros(A.nRings,numel(A.RunAnaObj.RunList));
    rateE = zeros(A.nRings,numel(A.RunAnaObj.RunList));
    cf    = zeros(A.nRings,numel(A.RunAnaObj.RunList));
    for i=1:A.nRings
        R           = A.MultiObj(i);
        count(i,:)  = R.SingleRunData.TBDIS_RM;
        sstime      = mean(R.SingleRunData.qUfrac_RM,1).*R.SingleRunData.TimeSec;
        rate(i,:)   = count(i,:)./sstime;
        rateE(i,:)  = sqrt(count(i,:))./sstime;
        cf(i,:)     = R.RMRateErosCorrectionqUActivity;
        
        Crate{j,i}                =  (rate(i,:).*cf(i,:)./Activity{j}.*ActivityRW2 - RefPeriodRW2(i));
        
        CrateEquivalentmV{j,i}    =  -(rate(i,:).*cf(i,:)./mean(Activity{j}).*ActivityRW2 - RefPeriodRW2(i)) ./737.8 * 1e3 * 117 / numel(R.PixList);
        
        rateEquivalentmV_E{j,i}   =  (rateE(i,:) ./737.8 *1e3 * 117 / numel(R.PixList));

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

save('SamakKNM2_DriftInRW123PSR1234_mVperDay.mat','SlopeRW123PSR1234_mV','SlopeErrorRW123PSR1234_mV');
<<<<<<< HEAD

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
=======
>>>>>>> 45c7fe312536697cee14a9669aa673f57829ba52
