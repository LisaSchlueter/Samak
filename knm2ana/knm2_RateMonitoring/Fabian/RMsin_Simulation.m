% KNM2 - -300V Rate Monitor FPD
% Stacked-pixel Per Pseudo-Ring: Simulation
% Conversion of rate in Potential Fluctuation
% Last Modified: 10/04/2020
% T. Lasserre
%

%% Read All Data - Loop on Period
for j=1:3
    
    RunList   = ['KNM2_RW' num2str(j)];
    
    % Read Data
    RunAnaArg = {'RunList',RunList,...
        'DataType','Real',...
        'FSDFlag','Sibille0p5eV',...
        'ELossFlag','KatrinT2A20',...
        'AnaFlag','StackPixel',...
        'RingMerge','Full',...
        'NonPoissonScaleFactor',1,...
        'ROIFlag','Default'};
    if ~exist('MR','var') || numel(MR)<3
        MR{j}        = MultiRunAnalysis(RunAnaArg{:});
        A{j}         = RingAnalysis('RunAnaObj',MR{j},'RingList',1:4);
    end
    
    % Stacked Pixel Data for each patch
    count  = zeros(A{j}.nRings,numel(A{j}.RunAnaObj.RunList));
    rate   = zeros(A{j}.nRings,numel(A{j}.RunAnaObj.RunList));
    rateE  = zeros(A{j}.nRings,numel(A{j}.RunAnaObj.RunList));
    sstime = zeros(A{j}.nRings,numel(A{j}.RunAnaObj.RunList));
    cf     = zeros(A{j}.nRings,numel(A{j}.RunAnaObj.RunList));
    
    for i=1:A{j}.nRings
        R           = A{j}.MultiObj(i);
        count(i,:)  = R.SingleRunData.TBDIS_RM;
        sstime(i,:) = mean(R.SingleRunData.qUfrac_RM,1).*R.SingleRunData.TimeSec;
        rate(i,:)   = count(i,:)./sstime(i,:);
        rateE(i,:)  = sqrt(count(i,:))./sstime(i,:);
        cf(i,:)     = R.RMRateErosCorrectionqUActivity;
        count_norm{j,i}     = rate(i,:) .* mean(sstime(i,:));
        corrcount_norm{j,i} = count_norm{j,i}  .* cf(i,:) ;
        scantime{j,i}       = sstime(i,:);
        RWperiodtime{j}     = sum(R.SingleRunData.qUfrac_RM.*R.SingleRunData.TimeSec);
        psrCount2eVSec{i}   =  -6.3 .* numel(R.PixList); % knm2
        StartTimeSecond{j}  = 1;
        ScanTimeSecond{j}   = seconds(R.SingleRunData.StartTimeStamp-R.SingleRunData.StartTimeStamp(1)+R.SingleRunData.TimeSec);
        StopTimeSecond{j}  = max(ScanTimeSecond{j});
    end
end

%% Start Simulation
for rwp          = 1:1 % period
    for psr          = 1:1 % ring
        
        % Timing Profile
        refscantime  = round(mean(scantime{rwp,psr}));       % average scan time per ring/period
        DeltaT       = 1;      % bin time, second
        ElapsedTime  = RWperiodtime{rwp}; % total time
        TimeArray    = 1:DeltaT:ElapsedTime;
        nTimeBins    = numel(TimeArray);
        
        % Reference Count Rates
        refcount     = mean(corrcount_norm{rwp,psr}); % average count per ring/period
        refcountE    = sqrt(refcount);                % uncertainty per ring/period
        
        % Add SINUSOIDAL Plasma Potential Fluctuations
        sigmaP           = 0.05;   % Volt
        PeriodModulation = 400;   % seconds
        %PeriodModulation = 60000;   % seconds
        CountModulation  = (psrCount2eVSec{psr} * sigmaP * DeltaT) .* sin( 2 * pi / PeriodModulation .* TimeArray + 2 * pi * rand(1));
        
        % Counts Simulation Per DeltaT
        CountDeltaTnaked = refcount./refscantime*DeltaT;
        CountDeltaT      = CountDeltaTnaked ...                          % reference rate / deltaT
            + randn(1,nTimeBins).*sqrt(CountDeltaTnaked); % statistical fluctuations
        CountDeltaT      = CountDeltaT + CountModulation; % Plasma
        
        % Counts Simulation Per KNM2 bin Time
        RMbinT        = floor(ElapsedTime./refscantime); % Number of knm2 scans
        CountRM       = (sum((reshape(CountDeltaT(1:RMbinT*refscantime),refscantime,RMbinT)),1))';
        CountRMmod    = (sum((reshape(CountModulation(1:RMbinT*refscantime),refscantime,RMbinT)),1))';

        displayPlot = 'ON';
        switch displayPlot
            case 'ON'
                
                %
                PlotStyle = { 's','MarkerSize',4,'MarkerFaceColor',rgb('CadetBlue'),...%rgb('IndianRed'),...%
                    'LineWidth',1,'Color',rgb('DarkBlue')};
                f22 = figure('Renderer','opengl');
                set(f22, 'Units', 'normalized', 'Position', [0.1, 0.1, 1, 0.8]);
                subplot(3,3,[1 2])
                plot(TimeArray,CountDeltaT,PlotStyle{:},'MarkerSize',2)
                hold on
                plot(TimeArray,CountDeltaTnaked.*ones(1,numel(TimeArray)),PlotStyle{:},'MarkerSize',2,'Color',rgb('Black'));
                hold off
                xlim([min(TimeArray) max(TimeArray)]);
                xlabel(sprintf('Time in %0.f sec bins [sec]',DeltaT));
                ylabel(sprintf('Counts / %0.f s',DeltaT));
                myMainTitle = sprintf('KATRIN - RW%0.f - FPD-PSR%0.f @E_0-300eV - Simulation',rwp,psr);
                title(myMainTitle);
                PrettyFigureFormat
                subplot(3,3,[4 5])
                p2=plot(TimeArray,CountModulation,'LineWidth',3,'Color','Black');
                leg=legend(p2,sprintf('Period %0.f seconds \n Amplitude %0.f mV',...
                    PeriodModulation,sigmaP*1000),'location','southeast');
                %leg.Color = 'none'; legend boxoff; 
                leg.Color='White';
                xlim([min(TimeArray) max(TimeArray)]);
                xlabel(sprintf('Time in %0.f sec bins [sec]',DeltaT));
                ylabel(sprintf('Modulation / %0.f s',DeltaT));
                PrettyFigureFormat
                subplot(3,3,[7 8])
                e1=errorbar(1:1:numel(CountRM),CountRM',sqrt(CountRM'),PlotStyle{:});
                e1.CapSize = 0;
                xlim([1 numel(CountRM)]);
                xlabel(sprintf('Scan Number RW%0.f',rwp));
                ylabel(sprintf('Counts / %0.f s',refscantime));
                PrettyFigureFormat
                subplot(3,3,[3 6 9])
                h=histogram(CountRM,15,'Normalization','pdf',...
                    'FaceColor',rgb('DarkBlue'),'LineWidth',2,'FaceAlpha',0.6);
                xlabel(sprintf('counts in %.1f sec',refscantime));
                ylabel('Frequency');
                PrettyFigureFormat
                
                % Gaussian PDF & fit
                pdfG = @(b) 1/sqrt(2*pi*pdG.sigma.^2) .* exp(-(b-pdG.mu).^2/(2*pdG.sigma.^2));
                pdG = fitdist(CountRM,'Normal');
                
                % Poisson PDF & fit
                pdfP = @(b) 1/sqrt(2*pi*pdN.lambda) .* exp(-(b-pdN.lambda).^2/(2*pdN.lambda));
                pdN = fitdist(CountRM,'poisson');
                
                hold on
                b = linspace(min(CountRM),max(CountRM),100);
                g=plot(b,pdfG(b),'Color',rgb('IndianRed'),'LineWidth',3);
                p=plot(b,pdfP(b),'Color',rgb('GoldenRod'),'LineWidth',3);
                leg=legend([h g p],...
                    sprintf('%.0f subscans',RMbinT),...
                    sprintf('\\sigma_g=%.0f counts',...
                    (pdG.sigma)),sprintf('\\sigma_p=%.0f counts \n Non-Poisson \n Factor = %.1f',...
                    sqrt(pdN.lambda),pdG.sigma/sqrt(pdN.lambda)),...
                    'location','northwest');
                leg.Color = 'none'; legend boxoff;
                hold off
                PrettyFigureFormat
                set(gca,'FontSize',12);
        end
        
    end
end

%% Toy MC - MultiPeriod - MultiRing - MultiToy
progressbar('Toy MC All RW All PSR');
counterRWPSR = 0;
for rwp          = 1:3 % period
    for psr          = 1:4 % ring
        % Toy MC of Many KNM2
        % progressbar('');
        Ntoy       = 500; NPnaked = zeros(1,Ntoy); NPplasma = zeros(1,Ntoy);
        
        Avector    = 0.0:0.025:0.5;
        Pvector    = decade(1,7,1);
        plasmaMap  = zeros(numel(Avector),numel(Pvector));
        
        counter  = 0;
        counterA = 0;
        for toyA = Avector
            counterA=counterA+1;
            counterP = 0;
            for toyP = Pvector
                counterP=counterP+1;
                counter = counter+1;
                for n=1:1:Ntoy
                    nakedC       = CountDeltaTnaked + randn(1,nTimeBins).*sqrt(CountDeltaTnaked);
                    Modulation   = (psrCount2eVSec{psr} * toyA * DeltaT) .* sin( 2 * pi / toyP .* TimeArray + 2 * pi * rand(1));
                    plasmaC      = nakedC + Modulation;
                    nakedCrm     = (sum((reshape(nakedC(1:RMbinT*refscantime),refscantime,RMbinT)),1))';
                    plasmaCrm    = (sum((reshape(plasmaC(1:RMbinT*refscantime),refscantime,RMbinT)),1))';
                    %pdGnaked(n)  = fitdist(nakedCrm,'Normal'); pdNnaked(n)  = fitdist(nakedCrm,'poisson');
                    pdGplasma(n) = fitdist(plasmaCrm,'Normal');pdNplasma(n) = fitdist(plasmaCrm,'poisson');
                    %NPnaked(n) = pdGnaked(n).sigma/sqrt(pdNnaked(n).lambda);
                    NPplasma(n)= pdGplasma(n).sigma/sqrt(pdNplasma(n).lambda);
                end
                plasmaMap(counterA,counterP)      = mean(NPplasma);
                plasmaMapSTD(counterA,counterP)   = std(NPplasma);
                
                % Final Results
                plasmaM{rwp,psr} = plasmaMap;
                plasmaS{rwp,psr} = plasmaMapSTD;
                
                % progressbar(counter/(numel(plasmaMap)));
            end
        end
        counterRWPSR=counterRWPSR+1;
        progressbar(counterRWPSR/numel(rwp)/numel(psr));
    end
end

%%
% figure(747)
% histogram(NPnaked)
% hold on
% histogram(NPplasma)
% hold off
% PrettyFigureFormat

%% Single Countour
myMainTitle = sprintf('KATRIN - RW%0.f - FPD-PSR%0.f @E_0-300eV - Simulation',rwp,psr);
fig1        = figure('Name',sprintf('KATRIN - %s - FPD @E_0-300eV - Simulation',RunList),...
    'NumberTitle','off','rend','painters','pos',[10 10 1000 1000]);
a           = annotation('textbox', [0 0.9 1 0.1], 'String', myMainTitle,'EdgeColor', ...
    'none','HorizontalAlignment', 'center','FontSize',22,'FontWeight','bold');
[X,Y] = meshgrid(Avector*1e3,log10(Pvector));
%[M,c] = contour(X,Y,plasmaM{3,2}',[1.1 1.2 1.3 1.4 1.5 2 2.5],'ShowText', 'on');
%c.LineWidth = 5;
clevels = [1.1:0.1:1.5 2:1:5];
hc = contourfcmap(X,Y,plasmaM{3,2}',clevels,flipud(winter(numel(clevels)-1)), ...
    'lo', rgb('White'), ...
    'hi', rgb('Navy'), ...
    'cbarloc', 'eastoutside', ...
    'method', 'calccontour');
hold on
line([0 500],[log10(MR{rwp}.RunData.TimeSec) log10(MR{rwp}.RunData.TimeSec)],'LineWidth',2,'LineStyle','--','Color',rgb('IndianRed'));
hold off
ylabel(hc.cb.ax, '\sigma_g / \sigma_p','FontSize',20);%\sigma_g / \sigma_p
xlabel('Amplitude [mV]'); xlim([0 500]);
ylabel('Log10(Period) [s]'); ylim([1 7]);
PrettyFigureFormat

%% All Contours - Results From Simulation
load('plasmaMapv2.mat');
myMainTitle = sprintf('KATRIN - All RW - All FPD-PSR @E_0-300eV - Simulation');
fig2        = figure('Name',sprintf('KATRIN - %s - FPD @E_0-300eV - Simulation',RunList),...
    'NumberTitle','off','rend','painters','pos',[10 10 1400 1000]);
a           = annotation('textbox', [0 0.9 1 0.1], 'String', myMainTitle,'EdgeColor', ...
    'none','HorizontalAlignment', 'center','FontSize',20,'FontWeight','bold');
counter=0;
for rwp=1:3
    for psr=1:4
        counter=counter+1;
        subplot(3,4,counter)
        clevels = [1.1:0.1:2.5];
        [X,Y] = meshgrid(Avector*1e3,log10(Pvector));
        hc = contourfcmap(X,Y,plasmaM{rwp,psr}',clevels,flipud(winter(numel(clevels)-1)), ...
            'lo', rgb('White'), ...
            'hi', rgb('Navy'), ...
            'cbarloc', 'eastoutside', ...
            'method', 'calccontour');
        hold on
        line([0 500],[log10(MR{rwp}.RunData.TimeSec) log10(MR{rwp}.RunData.TimeSec)],'LineWidth',2,'LineStyle','--','Color',rgb('IndianRed'));
        hold off
        ylabel(hc.cb.ax, '\sigma_g / \sigma_p','FontSize',12);%\sigma_g / \sigma_p
        PrettyFigureFormat
        xlabel(sprintf('RW%0.f - PSR%0.f - Amplitude [mV]',rwp,psr),'FontSize',12); xlim([0 500]);
        ylabel('Log10(Period) [s]','FontSize',12); ylim([1 7]);
    end
end

%% All Contours - Results From Simulation Compared to Data - Best Fit
load('plasmaMapv2.mat');
knm2dataE300V = [1.1 1.3 1.2 1.0; 1.1 1.1 1.0 1.0 ; 1.1 1.2 1.3 1.1];
myMainTitle = sprintf('KATRIN - All RW - All FPD-PSR @E_0-300eV - Data Best Fit Verus Simulation');
fig3        = figure('Name',sprintf('KATRIN - %s - FPD @E_0-300eV - Simulation',RunList),...
    'NumberTitle','off','rend','painters','pos',[10 10 1400 1000]);
a           = annotation('textbox', [0 0.9 1 0.1], 'String', myMainTitle,'EdgeColor', ...
    'none','HorizontalAlignment', 'center','FontSize',20,'FontWeight','bold');
counter=0;
for rwp=1:3
    for psr=1:4
        counter=counter+1;
        subplot(3,4,counter)
        clevels = [knm2dataE300V(rwp,psr)-0.005 knm2dataE300V(rwp,psr)+0.005];
        [X,Y] = meshgrid(Avector*1e3,log10(Pvector));
        hc = contourfcmap(X,Y,plasmaM{rwp,psr}',clevels,flipud(winter(numel(clevels)-1)), ...
            'lo', rgb('White'), ...
            'hi', rgb('White'), ...
            'cbarloc', 'eastoutside', ...
            'method', 'calccontour');
        hold on
        line([0 500],[log10(MR{rwp}.RunData.TimeSec) log10(MR{rwp}.RunData.TimeSec)],'LineWidth',2,'LineStyle','--','Color',rgb('IndianRed'));
        hold off
        ylabel(hc.cb.ax, '\sigma_g / \sigma_p','FontSize',12);%\sigma_g / \sigma_p
        PrettyFigureFormat
        xlabel(sprintf('RW%0.f - PSR%0.f - Amplitude [mV]',rwp,psr),'FontSize',12); xlim([0 500]);
        ylabel('Log10(Period) [s]','FontSize',12); ylim([1 7]);
    end
end

%% All Contours - Results From Simulation Compared to Data - Exclusion 2sigma
load('plasmaMapv2.mat');
myMainTitle = sprintf('KATRIN - All RW - All FPD-PSR @E_0-300eV - Data Best Fit Verus Simulation');
fig3        = figure('Name',sprintf('KATRIN - %s - FPD @E_0-300eV - Simulation',RunList),...
    'NumberTitle','off','rend','painters','pos',[10 10 1400 1000]);
a           = annotation('textbox', [0 0.9 1 0.1], 'String', myMainTitle,'EdgeColor', ...
    'none','HorizontalAlignment', 'center','FontSize',20,'FontWeight','bold');
counter=0;
for rwp=1:3
    for psr=1:4
        counter=counter+1;
        subplot(3,4,counter)
        clevels = [knm2dataE300V(rwp,psr)-0.005 knm2dataE300V(rwp,psr)+0.005];
[X,Y] = meshgrid(Avector*1e3,log10(Pvector));
hc = contourfcmap(X,Y,plasmaM{rwp,psr}'-2*plasmaS{rwp,psr}',clevels,flipud(winter(numel(clevels)-1)), ...
     'lo', rgb('SpringGreen'), ...
     'hi', rgb('LightCoral'), ...
     'cbarloc', 'eastoutside', ...
     'method', 'calccontour');
hold on
line([0 500],[log10(MR{rwp}.RunData.TimeSec) log10(MR{rwp}.RunData.TimeSec)],'LineWidth',2,'LineStyle','--');
hold off
ylabel(hc.cb.ax, '\sigma_g / \sigma_p','FontSize',12);%\sigma_g / \sigma_p
PrettyFigureFormat
xlabel(sprintf('RW%0.f - PSR%0.f - Amplitude [mV]',rwp,psr),'FontSize',12); xlim([0 500]);
ylabel('Log10(Period) [s]','FontSize',12); ylim([1 7]);
    end
end