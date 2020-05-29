Mode = 'Ringwise';
if strcmp(Mode,'Ringwise')
    load('SamakKNM2_CratemVequivalentInRW123PSR1234.mat');
    load('SamakKNM2_CratemVequivalentSlopesAndFitsInRW123PSR1234.mat');
    load('SamakKNM2_ShiftDriftInRW123PSR1234_mVperDay.mat');
    for j=1:3
        pause(1);
        RunList   = ['KNM2_RW' num2str(j)];
        %% Rate Evolution --> mV equivalent
        myMainTitle = sprintf('KATRIN - KNM2 RW%s - FPD Rate Evolution 300eV below Endpoint',num2str(j));
        maintitle   = myMainTitle;
        savefile1    = sprintf('plots/KNM2_RM300_EffectivePotentialFluctuation_RW%.0f_PseudoRings_1.png',j);
        fig1      = figure('Name',sprintf('KATRIN - %s Scanwise Background',RunList),...
            'NumberTitle','off','rend','painters','pos',[10 10 1400 1000]);
        a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
        a.FontSize=24;a.FontWeight='bold';

        for i=1:4

            subplot(4,1,i)
            hi=plot(StartTimeStampDays{j},CrateEquivalentmV{j,i},'s--','Color',rgb('IndianRed'),'LineWidth',1,'MarkerSize',12,'markerfacecolor',rgb('IndianRed'));

            hold on
            hc = plot(StartTimeStampDays{j},repmat(mean(CrateEquivalentmV{j,i}),1,numel(StartTimeStampDays{j})),'--','Color',rgb('Black'),'LineWidth',2);
            hold off

            hold on
            hl = plot(StartTimeStampDays{j},LinRW123PSR1234{j,i},'-','Color',rgb('DarkBlue'),'LineWidth',2);
            hold off

            hold on
            hq = plot(StartTimeStampDays{j},...
                Poly3RW123PSR1234{j,i},...
                '-','Color',rgb('DarkGreen'),'LineWidth',2);
            hold off

            ylabel('\Delta U (meV)');
            xlabel('Days since last RW setting');
            leg=legend([hi hc hl hq],...
                sprintf('Pseudo-Ring %0.f - std: %.2g meV ',i,std(CrateEquivalentmV{j,i})),...
                sprintf('constant - p-val=%.2g',pValCons{j,i}),...
                sprintf('linear: %.1f+-%.1f mV/day - p-val=%.2g',SlopeRW123PSR1234_mV{j,i},SlopeErrorRW123PSR1234_mV{j,i},pValLin{j,i}),...
                sprintf('poly3 - p-val=%.2g',pValPoly3{j,i}),...
                'location','eastoutside','FontSize',12);
            leg.Color = 'none'; legend boxoff;
            PrettyFigureFormat
        end
        
%         delete(hi);
%         delete(hc);
%         delete(hl);
%         delete(hq);

        %% Rate Evolution
        myMainTitle = sprintf('KATRIN - KNM2 RW%s - FPD Rate Evolution 300eV below Endpoint',num2str(j));
        maintitle   = myMainTitle;
        savefile2    = sprintf('plots/KNM2_RM300_EffectivePotentialFluctuation_RW%.0f_PseudoRings_2.png',j);
        fig2      = figure('Name',sprintf('KATRIN - %s Scanwise Background',RunList),...
            'NumberTitle','off','rend','painters','pos',[10 10 1400 1000]);
        a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
        a.FontSize=24;a.FontWeight='bold';

        for i=1:4

            subplot(4,1,i)
            hnc=plot(OverallStartTimeStamp{j},rateRaw{j,i},'s--','Color',rgb('Green'),'LineWidth',1,'MarkerSize',8,'markerfacecolor',rgb('DarkGreen'));
            hold on
            hc=plot(OverallStartTimeStamp{j},rate{j,i},'s--','Color',rgb('Red'),'LineWidth',1,'MarkerSize',12,'markerfacecolor',rgb('IndianRed'));
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
    pause(1);
      myMainTitle = sprintf('KATRIN KNM2 - FPD Rate 300eV below Endpoint - Average WGTS Potential Drift');
        maintitle   = myMainTitle;
        savefile3    = sprintf('plots/KNM2_RM300_EffectivePotentialFluctuation_RW%.0f_PseudoRings_3.png',123);
        fig3      = figure('Name',sprintf('KATRIN - %s Scanwise Background',RunList),...
            'NumberTitle','off','rend','painters','pos',[10 10 1400 1000]);
        a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
        a.FontSize=24;a.FontWeight='bold';


        % ABSOLUTE SHIFTS CORRECTION

    for i=1:4
        TimeRW123        = [OverallStartTimeStamp{1} OverallStartTimeStamp{2} OverallStartTimeStamp{3}];
        RatemVRW123      = [CrateEquivalentmV{1,i} CrateEquivalentmV{2,i} CrateEquivalentmV{3,i}];
        RateErrormVRW123 = [rateEquivalentmV_E{1,i} rateEquivalentmV_E{2,i} rateEquivalentmV_E{3,i}];
        FitsRW123        = [LinRW123PSR1234{1,i} LinRW123PSR1234{2,i} LinRW123PSR1234{3,i}];

        subplot(4,1,i)
        H1=plot(TimeRW123,RatemVRW123,'s--','Color',rgb('Green'),'LineWidth',1,'MarkerSize',8,'markerfacecolor',rgb('IndianRed'));

        % Fit: Linear
        fitType = fittype('a*x + b'); % The equation for your fit goes here
        [fLin,gofLin] = fit((days(TimeRW123-TimeRW123(1)))',RatemVRW123',fitType,...
            'Weights',(1./RateErrormVRW123').^2,...
            'StartPoint',[0 RatemVRW123(1)],...
            'Robust','Bisquare');
        ci = confint(fLin,0.68); uncertaintyLin = (ci(2,1)-ci(1,1))/2;
        hold on
        H2 = plot(OverallStartTimeStamp{1} ,LinRW123PSR1234{1,i},'-','Color',rgb('DarkBlue'),'LineWidth',2);
        H3 = plot(OverallStartTimeStamp{2} ,LinRW123PSR1234{2,i},'-','Color',rgb('DarkRed'),'LineWidth',2);
        H4 = plot(OverallStartTimeStamp{3} ,LinRW123PSR1234{3,i},'-','Color',rgb('DarkOrange'),'LineWidth',2);
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


%     %% Correction Difference between methods
%     % OverallStartTimeStamp
%     % CrateEquivalentmV
%     % rateEquivalentmV_E
% 
%       myMainTitle = sprintf('KATRIN KNM2 - Correction Difference');
%         maintitle   = myMainTitle;
%         savefile3    = sprintf('plots/KNM2_RM300_CorrectionDiff_RW%.0f_PseudoRings_3.png',123);
%         fig4      = figure('Name',sprintf('KATRIN - %s Correction Diff',RunList),...
%             'NumberTitle','off','rend','painters','pos',[10 10 1400 1000]);
%         a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
%         a.FontSize=24;a.FontWeight='bold';
% 
% 
%         % ABSOLUTE SHIFTS CORRECTION
% 
%     for i=1:4
%         TimeRW123        = [OverallStartTimeStamp{1} OverallStartTimeStamp{2} OverallStartTimeStamp{3}];
%         Rate1RW123       = [Crate1{1,i} Crate1{2,i} Crate1{3,i}];
%         Rate2RW123       = [Crate2{1,i} Crate2{2,i} Crate2{3,i}];
%         percentages      = std(Rate2RW123-Rate1RW123)./mean(rateRaw(:,i));
% 
%         subplot(4,1,i)
%         H1=plot(TimeRW123,(Rate2RW123-Rate1RW123),'s--','Color',rgb('Green'),'LineWidth',1,'MarkerSize',8,'markerfacecolor',rgb('IndianRed'));
% 
%         ylabel('(\Delta r)/r (%)');
%         ylim([-inf inf]);
%             xlabel('Days');
%             leg=legend([H1],...
%                 sprintf('Pseudo-Ring %0.f',i),'location','southeast','FontSize',12);
% 
%             leg.Color = 'none'; legend boxoff;
%             PrettyFigureFormat
%     end
    
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
    ht=title(sprintf('Mean Potential Per Period Per Pseudo-ring'));
    ylabel('RW Period'); ylim([0.5 3.5]) ; yticks([1 2 3 4]); %ylabel({'RW1';'RW2';'RW3'});
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
    ht=title(sprintf('Slope Per Period Per Pseudo-ring'));
    ylabel('RW Period'); ylim([0.5 3.5]) ; yticks([1 2 3 4]); %ylabel({'RW1';'RW2';'RW3'});
    xlabel('Pseudo-Ring'); xlim([0.5 4.5]); xticks([1 2 3 4]);%xlabel({'PSR1','PSR2','PSR3','PSR4'});
    zlabel('mV-equivalent / day'); %zlim([0 10]);
    PrettyFigureFormat
    ht.FontSize=16;
end