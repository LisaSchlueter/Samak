% ------------------------------------------------------------
% KNM2 rate monitor point (-300eV) analysis
% calculate rate - qU dependence with MC simulations
% ------------------------------------------------------------

% You need to have Run GetSamakPath in your Working Directory Before
savedir = [getenv('SamakPath'),'knm2ana/knm2_RateMonitoring/results'];
MakeDir(savedir);

% Fake KNM2 Data - RW2 Confirguration
% All KNM2 RW2 inputs are in ref_FakeRun_KNM2_CD84_RateMonitor
% (fake run with 20 data points +-1eV around E0-300eV)
savename = [savedir,'knm2_RateMonitoring_qUSlope.mat'];
RecomputeFlag = 'ON';
if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename);
else
    InitFile =  @ref_FakeRun_KNM2_CD84_RateMonitor; 
    CommonArg = {'FakeInitFile',InitFile,...
        'RunNr',1,...% has no meaning
        'DataType','Fake',...
        'FSDFlag','BlindingKNM2',...
        'ELossFlag','KatrinT2',...
        'chi2','chi2Stat',...
        'RingMerge','Full',...
        'minuitOpt','min;migrad',...
        'NonPoissonScaleFactor',1,...
        'AnaFlag','StackPixel',...
        'fixPar','mNu E0 Bkg Norm',...
        'RingList',1:12};
    % set up RunAnalysis object:
    A = RunAnalysis(CommonArg{:});
    % set ring segmentation
    R = RingAnalysis('RunAnaObj',A,'RingList',1:4);

    % get -300 eV rate / error
    qUIndex = find((A.RunData.qU-18574)<-95); % all Data points around rate monitor point
    qU      = A.RunData.qU(qUIndex);
    Rate    = A.RunData.TBDIS(qUIndex)./(A.RunData.qUfrac().*A.RunData.TimeSec);
    RateErr = sqrt(Rate);
     
    % linear fit of rate around E0-300V
    [par, err, chi2min,dof] = linFit(qU-18574,Rate,RateErr);
    
    save(savename,'par','err','chi2min','dof','qUIndex','qU','Rate','RateErr','CommonArg');
end

%% Load KNM2 Period RW2 REAL data for comparison
%  Normalization may be slightly different
    Knm2AnaArg = {'RunList','KNM2_RW2','DataType','Real',...
        'FSDFlag','BlindingKNM2','ELossFlag','KatrinT2',...
        'AnaFlag','StackPixel','RingMerge','Full','AnaFlag','StackPixel'};
    D  = MultiRunAnalysis(Knm2AnaArg{:});
    DR = RingAnalysis('RunAnaObj',D,'RingList',1:4);
    MeanRateDataUniform    = mean(D.SingleRunData.TBDIS_RM)./(mean(D.SingleRunData.qUfrac_RM).*mean(D.SingleRunData.TimeSec));
    qUDataUniform          = mean(D.SingleRunData.qU_RM);
        
%% Sanity plot: Uniform - All Pixels
f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
e1 = errorbar(qU-18574,Rate,RateErr,'-x','LineWidth',2);
hold on;
x = linspace(min(qU),max(qU),100)-18574;
y = par(1).*x+par(2);
p1 = plot(x,y,'LineWidth',2);
ud=plot(qUDataUniform-18574,MeanRateDataUniform,'rs','MarkerSize',20);
hold off;
PrettyFigureFormat('FontSize',22);
xlabel('Retarding energy - 18574 (eV)');
ylabel('Rate (cps)');
leg = legend([e1,p1,ud],sprintf('MC simulation'),...
    sprintf('Linear fit slope = %.2f \\pm %.2f cps/mV',par(1)/1e3,err(1)/1e3),...
    'KNM2 RW2 Data');
leg.EdgeColor = rgb('Silver');

%% Sanity plot: Pseudo-Ring Wise
    clear par err;
    myMainTitle = sprintf('KATRIN');
    maintitle   = myMainTitle;
    savefile    = sprintf('plots/tmp.png');
    fig      = figure('Name',sprintf('KATRIN'),...
        'NumberTitle','off','rend','painters','pos',[10 10 1400 1000]);
    a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
    a.FontSize=24;a.FontWeight='bold';
    
    qUIndex      = find((R.MultiObj(1).RunData.qU-18574)<-95); % all Data points around rate monitor point
    qU      = zeros(A.nRings,numel(qUIndex));
    Rate    = zeros(A.nRings,numel(qUIndex));
    RateErr = zeros(A.nRings,numel(qUIndex));
    
    for i=1:A.nRings
        
        CurrentRing = R.MultiObj(i);
         
        sp(i) = subplot(1,A.nRings,i)
        
        % Data Per Ring
        CurrentRingData = DR.MultiObj(i);
        MeanRateData(i) = mean(CurrentRingData.SingleRunData.TBDIS_RM)./(mean(CurrentRingData.SingleRunData.qUfrac_RM).*mean(CurrentRingData.SingleRunData.TimeSec));
        qUData(i)       = mean(CurrentRingData.SingleRunData.qU_RM);
    
        % get -300 eV points
        qU(i,:)      = CurrentRing.RunData.qU(qUIndex)';
        Rate(i,:)    = CurrentRing.RunData.TBDIS(qUIndex)'./(CurrentRing.RunData.qUfrac().*CurrentRing.RunData.TimeSec);
        RateErr(i,:) = sqrt(Rate(i,:));
        
        % linear fit
        [tmppar, tmperr, chi2min,dof] = linFit(qU(i,:)'-18574,Rate(i,:)',RateErr(i,:)');
        parRing(i,:)=tmppar;
        errRing(i,:)=tmperr;
        
        % Plot
        e1 = errorbar(qU(i,:)-18574,Rate(i,:),RateErr(i,:).*0,'-x','LineWidth',5);
        hold on;
        x = linspace(min(qU(i,:)),max(qU(i,:)),100)-18574;
        y = parRing(i,1).*x + parRing(i,2);
        p1 = plot(x,y,'LineWidth',5);
        rd = plot(qUData(i)-18574,MeanRateData(i),'rs','MarkerSize',20,'LineWidth',4);
        hold off;
        PrettyFigureFormat('FontSize',22);
        xlabel('Retarding energy - 18574 (eV)');
        ylabel('Rate (cps)');
        leg = legend([e1,p1,rd],sprintf('MC simulation'),...
            sprintf('Linear fit slope = %.2f \\pm %.2f cps/mV',parRing(i,1)/1e3,errRing(i,1)/1e3),...
            sprintf('KNM2 RW2 Data - PSR %0.f',i));
        leg.EdgeColor = rgb('Silver');
        leg.FontSize  = 8;
    end
    linkaxes([sp(1),sp(2),sp(3),sp(4)],'x');
    xlim([-300 -299.5])

    
 %% Renormalization PSR1
MeanRateData(1);
MeanRateSim = parRing(1,1).*(qUData(1)-18574) + parRing(1,2);
CorrectionFactorPSR1 = MeanRateSim/MeanRateData(1);

    myMainTitle = sprintf('KATRIN');
    maintitle   = myMainTitle;
    savefile    = sprintf('plots/tmp.png');
    fig      = figure('Name',sprintf('KATRIN'),...
        'NumberTitle','off','rend','painters','pos',[10 10 1400 1000]);
    a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
    a.FontSize=24;a.FontWeight='bold';
    
    qUIndex      = find((R.MultiObj(1).RunData.qU-18574)<-95); % all Data points around rate monitor point
    qU      = zeros(A.nRings,numel(qUIndex));
    Rate    = zeros(A.nRings,numel(qUIndex));
    RateErr = zeros(A.nRings,numel(qUIndex));
    
    for i=1:A.nRings
        
        CurrentRing = R.MultiObj(i);
         
        sp(i) = subplot(1,A.nRings,i)
        
        % Data Per Ring
        CurrentRingData = DR.MultiObj(i);
        MeanRateData(i) = mean(CurrentRingData.SingleRunData.TBDIS_RM)./(mean(CurrentRingData.SingleRunData.qUfrac_RM).*mean(CurrentRingData.SingleRunData.TimeSec));
        qUData(i)       = mean(CurrentRingData.SingleRunData.qU_RM);
    
        % get -300 eV points
        qU(i,:)      = CurrentRing.RunData.qU(qUIndex)';
        Rate(i,:)    = CurrentRing.RunData.TBDIS(qUIndex)'./(CurrentRing.RunData.qUfrac().*CurrentRing.RunData.TimeSec);
        RateErr(i,:) = sqrt(Rate(i,:));
        
        % Plot
        e1 = errorbar(qU(i,:)-18574,Rate(i,:),RateErr(i,:).*0,'-x','LineWidth',5);
        hold on;
        x = linspace(min(qU(i,:)),max(qU(i,:)),100)-18574;
        y = parRing(i,1).*x + parRing(i,2);
        p1 = plot(x,y,'LineWidth',5);
        rd = plot(qUData(i)-18574,MeanRateData(i).*CorrectionFactorPSR1,'rs','MarkerSize',20,'LineWidth',4);
        hold off;
        PrettyFigureFormat('FontSize',22);
        xlabel('Retarding energy - 18574 (eV)');
        ylabel('Rate (cps)');
        leg = legend([e1,p1,rd],sprintf('MC simulation'),...
            sprintf('Linear fit slope = %.2f \\pm %.2f cps/mV',parRing(i,1)/1e3,errRing(i,1)/1e3),...
            sprintf('KNM2 RW2 Data - PSR %0.f',i));
        leg.EdgeColor = rgb('Silver');
        leg.FontSize  = 8;
    end
    linkaxes([sp(1),sp(2),sp(3),sp(4)],'x');
    xlim([-301 -299])
    
    %% Calculate Rate Differences & Convert into Potential Shifts
    for i=1:4
        RateDifferences(i) = (parRing(i,1).*(qUData(i)-18574) + parRing(i,2)) - MeanRateData(i).*CorrectionFactorPSR1;
        mVDifferences(i)   = RateDifferences(i)./(parRing(i,1)/1e3);
    end
    
    