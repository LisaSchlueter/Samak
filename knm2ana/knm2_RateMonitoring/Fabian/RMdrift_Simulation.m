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
        % RWperiodtime{j}     = sum(R.SingleRunData.qUfrac_RM.*R.SingleRunData.TimeSec);
        psrCount2eVSec{i}   =  -6.3 .* numel(R.PixList); % knm2
        StartTimeSecond{j}  = 1;
        ScanTimeSecond{j}   = seconds(R.SingleRunData.StartTimeStamp-R.SingleRunData.StartTimeStamp(1));
        StopTimeSecond{j}   = max(ScanTimeSecond{j});
        RWperiodtime{j}     = StopTimeSecond{j} + R.SingleRunData.TimeSec(end);
    end
end

%% Start Simulation
for rwp          = 2:2 % period
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
        
        % Time Renormalization
        
        % Add SINUSOIDAL Plasma Potential Fluctuations
        VperDay          = 1*0.006;   % Volt, renormalized to Actual Period Duration
        CountDrift       = (psrCount2eVSec{psr} * VperDay/86400 * DeltaT)  .* TimeArray ;
        
        % Counts Simulation Per DeltaT
        CountDeltaTnaked = refcount./refscantime*DeltaT;
        CountDeltaT      = CountDeltaTnaked ...              % reference rate / deltaT
            + randn(1,nTimeBins).*sqrt(CountDeltaTnaked);    % statistical fluctuations
        CountDeltaT      = CountDeltaT + CountDrift;         % Plasma
        
        % Counts Simulation Per KNM2 bin Time
        RMbinT         = numel(ScanTimeSecond{rwp}); % Number of knm2 scans
        %CountRM       = (sum((reshape(CountDeltaT(1:RMbinT*refscantime),refscantime,RMbinT)),1))';
        %CountRMmod    = (sum((reshape(CountDrift(1:RMbinT*refscantime),refscantime,RMbinT)),1))';
        for scan = 1:1:RMbinT
        CountRM(scan)  = sum(CountDeltaT(1+ScanTimeSecond{rwp}(scan):ScanTimeSecond{rwp}(scan)+refscantime));
        CountRMmod(scan)  = sum(CountDrift(1+ScanTimeSecond{rwp}(scan):ScanTimeSecond{rwp}(scan)+refscantime));
        end
        
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
                p2=plot(TimeArray,CountDrift,'LineWidth',3,'Color','Black');
                leg=legend(p2,sprintf('Drift %0.f mV/day',VperDay*1e3),'location','southWest');
                leg.Color = 'none'; legend boxoff; 
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
                pdG = fitdist(CountRM','Normal');
                
                % Poisson PDF & fit
                pdfP = @(b) 1/sqrt(2*pi*pdN.lambda) .* exp(-(b-pdN.lambda).^2/(2*pdN.lambda));
                pdN = fitdist(CountRM','poisson');
                
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
progressbar('Toy Drift MC All RW All PSR');
counterRWPSR = 0;
for rwp          = 1:3 % period
    for psr          = 1:4 % ring
        % Toy MC of Many KNM2
        Ntoy       = 250; NPnaked = zeros(1,Ntoy); NPplasma = zeros(1,Ntoy);
        
        Dvector         = [0 1 2 3 4 5 6 7 8 9 10 15 20]*1e-3;
        plasmaDriftMap  = zeros(numel(Dvector),1);
        counterD = 0;
        
        CountDeltaTnaked = mean(corrcount_norm{rwp,psr})./refscantime*DeltaT;

        for toyD = Dvector
            counterD= counterD+1;
            for n=1:1:Ntoy
                
                nakedC       = CountDeltaTnaked + randn(1,numel(1:DeltaT:RWperiodtime{rwp})).*sqrt(CountDeltaTnaked);
                CountDrift   = (psrCount2eVSec{psr} * toyD /86400 * DeltaT)  .* (1:DeltaT:RWperiodtime{rwp});
                plasmaC      = nakedC + CountDrift;
                
                %nakedCrm     = (sum((reshape(nakedC(1:RMbinT*refscantime),refscantime,RMbinT)),1))';
                %plasmaCrm    = (sum((reshape(plasmaC(1:RMbinT*refscantime),refscantime,RMbinT)),1))';
                RMbinT        = numel(ScanTimeSecond{rwp}); % Number of knm2 scans
                nakedCrm      = zeros(1,RMbinT);
                plasmaCrm     = zeros(1,RMbinT);
                for scan = 1:1:RMbinT
                    nakedCrm(scan)  = sum(nakedC(1+ScanTimeSecond{rwp}(scan):ScanTimeSecond{rwp}(scan)+refscantime));
                    plasmaCrm(scan) = sum(plasmaC(1+ScanTimeSecond{rwp}(scan):ScanTimeSecond{rwp}(scan)+refscantime));
                end
                
                pdGplasma    = fitdist(plasmaCrm','Normal');pdNplasma  = fitdist(plasmaCrm','poisson');
                NPplasma(n)  = pdGplasma.sigma/sqrt(pdNplasma.lambda);
            end
            plasmaDriftMap(counterD)      = mean(NPplasma);
            plasmaDriftMapSTD(counterD)   = std(NPplasma);
            
            % Final Results
            plasmaDriftM{rwp,psr} = plasmaDriftMap;
            plasmaDriftS{rwp,psr} = plasmaDriftMapSTD;
            
        end
    counterRWPSR=counterRWPSR+1;
    progressbar(counterRWPSR/12);
    end
end

%% All RW / PSR - Results From Simulation
myMainTitle = sprintf('KATRIN - All RW - All FPD-PSR @E_0-300eV - Drift Simulation');
fig2        = figure('Name',sprintf('KATRIN - %s - FPD @E_0-300eV - Drift Simulation',RunList),...
    'NumberTitle','off','rend','painters','pos',[10 10 1400 1000]);
a           = annotation('textbox', [0 0.9 1 0.1], 'String', myMainTitle,'EdgeColor', ...
    'none','HorizontalAlignment', 'center','FontSize',20,'FontWeight','bold');
counter=0;
for rwp=1:3
    for psr=1:4
        counter=counter+1;
        subplot(3,4,counter)
        plot(Dvector*1e3,(smooth(plasmaDriftM{rwp,psr},'loess')),'LineWidth',3,'Color',rgb('DarkBlue'));
        PrettyFigureFormat
        ylabel('\sigma_g / \sigma_p','FontSize',12); 
        xlabel(sprintf('RW%0.f - PSR%0.f - Drift [mV/Day]',rwp,psr),'FontSize',12);
    end
end