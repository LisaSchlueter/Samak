% KNM2 - -300V Rate Monitor FPD
% Stakced-pixel Evolution
% Pseudo-ring Wise
% Calculate Mean Rate + Error 
% Per Pseudo Ring 1-2-3-4
% Per Perdiod 1 2 3
%
% Last Modified: 04/02/2020
% T. Lasserre
%

% ServerConfig
LocalPath=getenv('SamakPath');
LocalPath=[LocalPath '/knm2ana/knm2_RateMonitoring'];
MakeDir(LocalPath);
%

MeanRate     = zeros(3,5);
SEMRate      = zeros(3,5);
MeanEffU     = zeros(3,5);
SEMEffU      = zeros(3,5);

%% Period 2 Mean Rate
%  Period 2 is stable, therefore i use it to normalize reference rate, per
%  pseudo ring. Period 1 and 3 are compared to period 2
RefPeriodRW2 =    1.0e+04 * [1.576410807682920 1.174928552130216 0.951121987964973 1.341066577847976 1.342891737873775];
ActivityRW2  =    4.193018127856629e+17;

for j=1:3
    
    RunList   = ['KNM2_RW' num2str(j)];
    %% Read Data
    DataType  = 'Real';
    FSDFlag   = 'BlindingKNM2';
    ELossFlag = 'KatrinT2';
    AnaFlag   = 'StackPixel'; % uniform FPD
    RunAnaArg = {'RunList',RunList,'DataType',DataType,...
        'FSDFlag',FSDFlag,'ELossFlag',ELossFlag,'AnaFlag',AnaFlag,'RingMerge','Azi'};
    
    MR        = MultiRunAnalysis(RunAnaArg{:});
    A         = RingAnalysis('RunAnaObj',MR,'RingList',1:5);
    R         = A.MultiObj(1);
    
    %% Slow Control Data
    p1 =(R.SingleRunData.WGTS_MolFrac_TT'+0.5*R.SingleRunData.WGTS_MolFrac_HT'+0.5*R.SingleRunData.WGTS_MolFrac_DT')./mean((R.SingleRunData.WGTS_MolFrac_TT'+0.5*R.SingleRunData.WGTS_MolFrac_HT'+0.5*R.SingleRunData.WGTS_MolFrac_DT')).*R.SingleRunData.WGTS_CD_MolPerCm2'./mean(R.SingleRunData.WGTS_CD_MolPerCm2);
    p2 = mean(R.SingleRunData.qU_RM,1); p2=p2-mean(p2);
    
    Activity  =  (R.SingleRunData.WGTS_MolFrac_TT'+0.5*R.SingleRunData.WGTS_MolFrac_HT'+0.5*R.SingleRunData.WGTS_MolFrac_DT').*(R.SingleRunData.WGTS_CD_MolPerCm2');

    %% Time in days
    StartTimeStampDays = days(R.SingleRunData.StartTimeStamp-R.SingleRunData.StartTimeStamp(1));
    
    %% Stacked Pixel Data for each patch
    count       = zeros(A.nRings,numel(A.RunAnaObj.RunList));
    rate        = zeros(A.nRings,numel(A.RunAnaObj.RunList));
    cf          = zeros(A.nRings,numel(A.RunAnaObj.RunList));
    
    for i=1:A.nRings
        R             = A.MultiObj(i);
        
        % Data
        count(i,:)    = R.SingleRunData.TBDIS_RM;
        sstime        = mean(R.SingleRunData.qUfrac_RM,1).*R.SingleRunData.TimeSec;
        
        % Rates
        rate(i,:)     = count(i,:)./sstime;
        cf(i,:)       = R.RMRateErosCorrectionqUActivity;
        
        % Rates Mean & SME
        MeanRate(j,i) = mean(rate(i,:).*cf(i,:).*ActivityRW2./Activity');
        SEMRate(j,i)  = std(rate(i,:) .*cf(i,:));
        
        % Effective Potential Shift RefPeriodRW2
        MeanEffU(j,i) = - (MeanRate(j,i) - RefPeriodRW2(i)) ./737.8 * 1e3 * 117 / numel(R.PixList); % mV
        
       SEMEffU(j,i)  =   SEMRate(j,i) ./ 737.8 * 1e3 * 117 / numel(R.PixList); %mV

    end
end

%% Save Data for Samak Analysis
ShiftRW12AziPatch12345 =  MeanEffU(2,:) - MeanEffU(1,:);
ShiftRW23AziPatch12345 =  MeanEffU(3,:) - MeanEffU(2,:);
ShiftErrorRW12AziPatch12345 = sqrt(SEMEffU(1,:).^2 + SEMEffU(2,:).^2);
ShiftErrorRW23AziPatch12345 = sqrt(SEMEffU(2,:).^2 + SEMEffU(3,:).^2);
save('SamakKNM2_ShiftRW123AziPatch12345_mV.mat','ShiftRW12PSR1234AziPatch12345','ShiftRW23AziPatch12345','ShiftErrorRW12AziPatch12345','ShiftErrorRW23AziPatch12345');


%% Plot Optional
PlotFlag = 'ON';
switch PlotFlag
    case 'ON'
        
        % Rates
        myMainTitle = sprintf('KATRIN - FPD Rate Evolution 300eV below Endpoint');
        maintitle   = myMainTitle;
        savefile    = sprintf('plots/test.png');
        fig1      = figure('Name',sprintf('KATRIN - %s Scanwise Background',RunList),...
            'NumberTitle','off','rend','painters','pos',[10 10 1400 1000]);
        a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
        a.FontSize=24;a.FontWeight='bold';
        
        for i=1:A.nRings
            subplot(5,2,2*i-1)
            errorbar([1 2 3],MeanRate(:,i),SEMRate(:,i),'s','MarkerSize',5,'LineWidth',2,'MarkerFaceColor',rgb('DarkBlue'),'Color',rgb('DarkBlue'));
            xticks([1 2 3]);
            xticklabels({'RW1','RW2','RW3'});
            xlim([0.9 3.1]);
            ylabel(sprintf('PR %0.f (cps)',i));
            PrettyFigureFormat
            
            subplot(5,2,2*i)
            errorbar([1 2 3],MeanEffU(:,i),SEMEffU(:,i),'s','MarkerSize',5,'LineWidth',2,'MarkerFaceColor',rgb('IndianRed'),'Color',rgb('IndianRed'));
            xticks([1 2 3]);
            xticklabels({'RW1','RW2','RW3'});
            xlim([0.9 3.1]);
            ylabel(sprintf('PR %0.f (mV)',i));
            ylim([-200 +200]);
            PrettyFigureFormat
        end
        
end