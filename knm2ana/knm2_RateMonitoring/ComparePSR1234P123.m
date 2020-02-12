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

%% Period 2 Mean Rate
%  Period 2 is stable, therefore i use it to normalize reference rate, per
%  pseudo ring. Period 1 and 3 are compared to period 2
RefPeriodRW2 =    1.0e+04 * [1.5764 2.0199 1.8992 1.0569];
ActivityRW2  =    4.193018127856629e+17;

%% Comparison with endpoint data
format long
E090RW1file = importdata([LocalPath '/endpointdata/SingleRingFitResult_Full_KNM2_RW1_171runs_freeParE0BkgNorm_chi2Stat_90eVrange.mat']);
E090RW2file = importdata([LocalPath '/endpointdata/SingleRingFitResult_Full_KNM2_RW2_109runs_freeParE0BkgNorm_chi2Stat_90eVrange.mat']);
E090RW3file = importdata([LocalPath '/endpointdata/SingleRingFitResult_Full_KNM2_RW3_84runs_freeParE0BkgNorm14_chi2Stat_90eVrange.mat']);

for i=1:4
E090RW1PR(i)   = (E090RW1file(1).par(i,2) - E090RW1file(1).par(i,2))*1e3;%+ 18574;
E090RW2PR(i)   = (E090RW2file(1).par(i,2) - E090RW1file(1).par(i,2) )*1e3;%+ 18574;
E090RW3PR(i)   = (E090RW3file(1).par(i,2) - E090RW1file(1).par(i,2) )*1e3;%+ 18574;

E090RW1PRe(i)  = (E090RW1file(1).err(i,2) - E090RW1file(1).par(i,2) )*1e3;%+ 18574;
E090RW2PRe(i)  = (E090RW2file(1).err(i,2) - E090RW1file(1).par(i,2) )*1e3;%+ 18574;
E090RW3PRe(i)  = (E090RW3file(1).err(i,2) - E090RW1file(1).par(i,2) )*1e3;%+ 18574;

end

for i=1:4
E090(1,i)      = E090RW1PR(i);
E090(2,i)      = E090RW2PR(i);
E090(3,i)      = E090RW3PR(i);

E090e(1,i)     = E090RW1PR(i);
E090e(2,i)     = E090RW2PR(i);
E090e(3,i)     = E090RW3PR(i);
end
%%


MeanRate     = zeros(3,4);
SEMRate      = zeros(3,4);
MeanEffU     = zeros(3,4);
SEMEffU      = zeros(3,4);

for j=1:3
    
    RunList   = ['KNM2_RW' num2str(j)];
    %% Read Data
    DataType  = 'Real';
    FSDFlag   = 'BlindingKNM2';
    ELossFlag = 'KatrinT2';
    AnaFlag   = 'StackPixel'; % uniform FPD
    RunAnaArg = {'RunList',RunList,'DataType',DataType,...
        'FSDFlag',FSDFlag,'ELossFlag',ELossFlag,'AnaFlag',AnaFlag,'RingMerge','Full'};
    
    MR        = MultiRunAnalysis(RunAnaArg{:});
    A         = RingAnalysis('RunAnaObj',MR,'RingList',1:4);
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
        
        % Correction
        cf(i,:)       = R.RMRateErosCorrectionqUActivity;
        
        % Rates Mean & SME
        MeanRate(j,i) = mean(rate(i,:).*cf(i,:).*ActivityRW2./Activity');
        %SEMRate(j,i)  = std(rate(i,:) .*cf(i,:)) / sqrt(numel(rate(i,:)));
        SEMRate(j,i)  = std(rate(i,:) .*cf(i,:));
        
        % Effective Potential Shift RefPeriodRW2
        % MeanEffU(j,i) = - (MeanRate(j,i) - MeanRate(1,i)) ./737.8 * 1e3 * 117 / numel(R.PixList); % mV
        MeanEffU(j,i) = - (MeanRate(j,i) - RefPeriodRW2(i)) ./737.8 * 1e3 * 117 / numel(R.PixList); % mV
        
        SEMEffU(j,i)  =   SEMRate(j,i) ./ 737.8 * 1e3 * 117 / numel(R.PixList) / sqrt(numel(rate(i,:))); %mV
    end
end

%% Save Data for Samak Analysis
ShiftRW12PSR1234 =  MeanEffU(2,:) - MeanEffU(1,:);
ShiftRW23PSR1234 =  MeanEffU(3,:) - MeanEffU(2,:);
ShiftErrorRW12PSR1234 = sqrt(SEMEffU(1,:).^2 + SEMEffU(2,:).^2);
ShiftErrorRW23PSR1234 = sqrt(SEMEffU(2,:).^2 + SEMEffU(3,:).^2);
save('SamakKNM2_ShiftRW123PSR1234_mV.mat','ShiftRW12PSR1234','ShiftRW23PSR1234','ShiftErrorRW12PSR1234','ShiftErrorRW23PSR1234');

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
            subplot(4,2,2*i-1)
            errorbar([1 2 3],MeanRate(:,i),SEMRate(:,i),'s','MarkerSize',5,'LineWidth',2,'MarkerFaceColor',rgb('DarkBlue'),'Color',rgb('DarkBlue'));
            xticks([1 2 3]);
            xticklabels({'RW1','RW2','RW3'});
            xlim([0.9 3.1]);
            ylabel(sprintf('PR %0.f (cps)',i));
            PrettyFigureFormat
        end
        
        % Effective Potential Shift - 4 plots
        myMainTitle = sprintf('KATRIN - FPD Rate Evolution 300eV below Endpoint');
        maintitle   = myMainTitle;
        savefile    = sprintf('plots/test.png');
        fig2      = figure('Name',sprintf('KATRIN - %s Scanwise Background',RunList),...
            'NumberTitle','off','rend','painters','pos',[10 10 1400 1000]);
        a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
        a.FontSize=24;a.FontWeight='bold';
        
        for i=1:A.nRings
            subplot(4,1,i)
            errorbar([1 2 3],MeanEffU(:,i),SEMEffU(:,i),'s','MarkerSize',20,'LineWidth',2,'MarkerFaceColor',rgb('IndianRed'))

            
           % subplot(4,2,2*i)
           % errorbar([1 2 3],MeanEffU(:,i),SEMEffU(:,i),'s','MarkerSize',5,'LineWidth',2,'MarkerFaceColor',rgb('IndianRed'),'Color',rgb('IndianRed'));
            xticks([1 2 3]);
            xticklabels({'RW1','RW2','RW3'});
            xlim([0.9 3.1]);
            ylabel(sprintf('PR %0.f (mV)',i));
            ylim([-200 +200]);
            PrettyFigureFormat
        end
        
        
        % Effective Potential Shift - 1 Plot
        myMainTitle = sprintf('KATRIN - FPD Rate Evolution 300eV below Endpoint');
        maintitle   = myMainTitle;
        savefile    = sprintf('plots/test.png');
        fig3      = figure('Name',sprintf('KATRIN - %s Scanwise Background',RunList),...
            'NumberTitle','off','rend','painters','pos',[10 10 1400 1000]);
        a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
        a.FontSize=24;a.FontWeight='bold';
        
        for i=1:A.nRings
            hold on
            
            if i==2 || i==3
                MyWidth = 4;
            else
                MyWidth = 2;
            end

            h(i)=errorbar([1 2 3],MeanEffU(:,i),SEMEffU(:,i),'s-','MarkerSize',5,'LineWidth',MyWidth)

          %  h(i)=errorbar([1 2 3]+0.01*i,MeanEffU(:,i),SEMEffU(:,i),'s-','MarkerSize',5,'LineWidth',MyWidth)

            
            xticks([1 2 3]);
            xticklabels({'RW1','RW2','RW3'});
            xlim([0.9 3.1]);
            ylabel('mV');

            ylim([-60 +210]);
           % ylim([-200 +200]);

        end
        hrw1=line([-0.5 1.5],[-50 -50]+50,'LineStyle','--','LineWidth',2,'Color','Black');
        hrw2=line([1.5 2.5],[-8 -8]-42+50,'LineStyle','--','LineWidth',2,'Color','Black');
        hrw3=line([2.5 3.5],[+193 +193]-42+50,'LineStyle','--','LineWidth',2,'Color','Black');
        hold off
        legend([h(1),h(2),h(3),h(4),hrw1,hrw2,hrw3,],...
            'Pseudo-Ring 1','Pseudo-Ring 2','Pseudo-Ring 3','Pseudo-Ring 4',...
            'RW1 Setting','RW2 Setting','RW3 Setting','Location','NorthWest');
        PrettyFigureFormat
        
         % Effective Potential Shift - 1 Plot - Comparison with Endpoints
        myMainTitle = sprintf('KATRIN - FPD Rate Evolution 300eV below Endpoint');
        maintitle   = myMainTitle;
        savefile    = sprintf('plots/test.png');
        fig5      = figure('Name',sprintf('KATRIN - %s Scanwise Background',RunList),...
            'NumberTitle','off','rend','painters','pos',[10 10 1400 1000]);
        a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
        a.FontSize=24;a.FontWeight='bold';
        
        for i=1:A.nRings
            hold on
            
            if i==2 || i==3
                MyWidth = 2;
            else
                MyWidth = 2;
            end
            h(i)=errorbar([1 2 3]+0.01*i,E090([1 2 3],i),E090e([1 2 3],i),'s-','MarkerSize',5,'LineWidth',MyWidth);
            
            xticks([1 2 3]);
            xticklabels({'RW1','RW2','RW3'});
            xlim([0.9 3.1]);
            ylabel('Relative Endpoint Shift wrt RW1(meV)');
            %ylim([-200 +100]);
        end
        hrw1=line([-0.5 1.5],-([-50 -50]+50),'LineStyle','--','LineWidth',2,'Color','Black');
        hrw2=line([1.5 2.5],-([-8 -8]-42+50),'LineStyle','--','LineWidth',2,'Color','Black');
        hrw3=line([2.5 3.5],-([+193 +193]-42+50),'LineStyle','--','LineWidth',2,'Color','Black');
        hold off
        legend([h(1),h(2),h(3),h(4),hrw1,hrw2,hrw3,],...
            'Pseudo-Ring 1','Pseudo-Ring 2','Pseudo-Ring 3','Pseudo-Ring 4',...
            'RW1 Setting','RW2 Setting','RW3 Setting','Location','NorthWest');
        PrettyFigureFormat
        
                % Effective Potential Shift - 1 Plot - Comparison with Endpoints
        myMainTitle = sprintf('KATRIN - FPD Rate Evolution 300eV below Endpoint');
        maintitle   = myMainTitle;
        savefile    = sprintf('plots/test.png');
        fig4      = figure('Name',sprintf('KATRIN - %s Scanwise Background',RunList),...
            'NumberTitle','off','rend','painters','pos',[10 10 1400 1000]);
        a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
        a.FontSize=24;a.FontWeight='bold';
        
        for i=1:A.nRings
            hold on
            
            if i==2 || i==3
                MyWidth = 2;
            else
                MyWidth = 2;
            end
            h(i)=errorbar([1 2 3]+0.01*i,E090([1 2 3],i)+MeanEffU([1 2 3],i),E090e([1 2 3],i),'s-','MarkerSize',5,'LineWidth',MyWidth);
            
            xticks([1 2 3]);
            xticklabels({'RW1','RW2','RW3'});
            xlim([0.9 3.1]);
            ylabel('\Delta U_{eff} + Relative Endpoint Shift (meV)');
            %ylim([-200 +100]);
        end
       
        hold off
        legend([h(1),h(2),h(3),h(4)],...
            'Pseudo-Ring 1','Pseudo-Ring 2','Pseudo-Ring 3','Pseudo-Ring 4',...
            'Location','NorthWest');
        PrettyFigureFormat
end

