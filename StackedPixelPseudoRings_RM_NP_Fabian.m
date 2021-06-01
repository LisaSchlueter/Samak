QAplots = 'OFF';
HVDrift = 'OFF';
Detrend = 'OFF';
NBoot   = 1000;
ROIstr  = 'Default';
Corr    = 'Fabian';
Mode    = 'Periodwise';

if strcmp(Mode,'Ringwise')
    %% Rate Evolution --> mV equivalent
    myMainTitle = sprintf('KATRIN - KNM2 - FPD @E_0-300eV - Non-Poissonian Components');
    maintitle   = myMainTitle;
    fig1      = figure('Name',sprintf('KATRIN - KNM2 Scanwise Background'),...
        'NumberTitle','off','rend','painters','pos',[10 10 1400 1000]);
    a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
    a.FontSize=24;a.FontWeight='bold';

    savefile1    = sprintf('plots/KNM2_RM300_EffectivePotentialFluctuation.0f_PseudoRings_NP.png');


    % Loop on Period
    for j=1:3


        RunList   = ['KNM2_RW' num2str(j)];

        % Read Data
        DataType  = 'Real';
        RunAnaArg = {'RunList',RunList,'DataType',DataType,...
            'FSDFlag','BlindingKNM2','ELossFlag','KatrinT2',...
            'AnaFlag','StackPixel','RingMerge','Full','NonPoissonScaleFactor',1};
        R        = MultiRunAnalysis(RunAnaArg{:});
        range = 40;               % fit range in eV below endpoint        
        R.exclDataStart = R.GetexclDataStart(range); % find correct data, where to cut spectrum
        A         = RingAnalysis('RunAnaObj',R,'RingList',1:4);
        R         = A.MultiObj(1);
        R.ROIFlag=ROIstr; R.SetROI;

        % Slow Control Data
        p1 = (R.SingleRunData.WGTS_MolFrac_TT'+0.5*R.SingleRunData.WGTS_MolFrac_HT'+0.5*R.SingleRunData.WGTS_MolFrac_DT')./mean((R.SingleRunData.WGTS_MolFrac_TT'+0.5*R.SingleRunData.WGTS_MolFrac_HT'+0.5*R.SingleRunData.WGTS_MolFrac_DT')).*R.SingleRunData.WGTS_CD_MolPerCm2'./mean(R.SingleRunData.WGTS_CD_MolPerCm2');
        p2 = mean(R.SingleRunData.qU_RM,1); p2=p2-mean(p2);

        %% HV Drift Correction
        FirstDayPeriod1 = datetime('02-Oct-2019 14:52:19');
        TimeLineDaysFirstDayPeriod1 = days(R.SingleRunData.StartTimeStamp-FirstDayPeriod1);
        if strcmp(HVDrift,'ON')
            HVdriftPerPixel = 1.5*(TimeLineDaysFirstDayPeriod1) * 6.3e-3;
        else
            HVdriftPerPixel = 0;
        end

        %% Stacked Pixel Data for each patch
        count = zeros(A.nRings,numel(A.RunAnaObj.RunList));
        rate  = zeros(A.nRings,numel(A.RunAnaObj.RunList));
        rateE = zeros(A.nRings,numel(A.RunAnaObj.RunList));
        sstime    = zeros(A.nRings,numel(A.RunAnaObj.RunList));
        cf    = zeros(A.nRings,numel(A.RunAnaObj.RunList));
        for i=1:A.nRings
              R           = A.MultiObj(i);
              R.ROIFlag=ROIstr; R.SetROI;
              if strcmp(Corr,'Fabian')
                R.RMCorrection('QAplots',QAplots);
              end
              count  = R.SingleRunData.TBDIS_RM;
              sstime = mean(R.SingleRunData.qUfrac_RM,1).*R.SingleRunData.TimeSec;
              rate   = count./sstime;
              if strcmp(Detrend,'ON')
                  rate = detrend(rate)+mean(rate);
              end
              if strcmp(Corr,'Fabian')
                corrcount_norm{j,i} = rate .* mean(sstime);
              elseif strcmp(Corr,'Thierry')
                  cf(i,:)     = R.RMRateErosCorrectionqUActivity;
                corrcount_norm{j,i} = rate.*cf(i,:) .* mean(sstime);
              end
              
            if j==1
              subplot(3,4,i+j-1);
            elseif j==2
              subplot(3,4,4+i);
            elseif j==3
              subplot(3,4,8+i);
            end

            %% Stacked Pixel: Non-Poissonian Component?
            lambda=2;
            sigma=1;
            mu=0.1;
            par4=1;
            x=-2:0.0001:2;
            p = [lambda,sigma,mu,par4,x];
            [N,edges] = histcounts(corrcount_norm{j,i},15,'Normalization','probability');
            means = zeros(1,numel(N));
            for n=1:(numel(edges)-1)
                means(n) = edges(n)+(edges(n+1)-edges(n))./2;
            end
            means = means - mean(edges);
            Poisson = @(lambda,x) lambda.^x.*exp(1).^(-lambda)./gamma(x+1);
            Gauss   = @(sigma,mu,x) exp(-0.5.*((x-mu)./sigma).^2)./(sigma.*sqrt(2.*pi));
            Convolution = @(lambda,sigma,mu,par4,x) par4.*conv(Poisson(lambda,x),Gauss(sigma,mu,x),'same');
            %Conv = fit(means',N',stairs(Convolution),'StartPoint',p,'lower',[1e5,0,1e5,0],'upper',[1e8,1e2,1e8,10]);
            %plot(Conv);

            pdG = fitdist(corrcount_norm{j,i}','Normal');
            pdN = fitdist(corrcount_norm{j,i}','poisson');
            mydata=corrcount_norm{j,i}(abs(corrcount_norm{j,i}-pdG.mu)<(3*pdG.sigma))';
            pdG = fitdist(mydata,'Normal');
            pdN = fitdist(mydata,'poisson');
            Broadening(j,i) = sqrt(pdG.sigma^2-pdN.lambda)/(mean(sstime)*737.8) * 1e3 * 117 / numel(R.PixList);
            NP(j,i)         = pdG.sigma/sqrt(pdN.lambda);

            %% Error calculation: Bootstrap
            mydata_rand=zeros(numel(mydata),1);
            Broadening_rand=zeros(NBoot,1);
            NP_rand=zeros(NBoot,1);
            for l=1:NBoot
                  for k=1:numel(mydata)
                      mydata_rand(k) = mydata(randi(numel(mydata)));
                  end
                  pdG_rand = fitdist(mydata_rand,'Normal');
                  pdN_rand = fitdist(mydata_rand,'poisson');
                  Broadening_rand(l) = sqrt(abs(pdG_rand.sigma^2-pdN_rand.lambda))/(mean(sstime)*737.8) * 1e3 * 117 / numel(R.PixList);
                  NP_rand(l)         = pdG_rand.sigma/sqrt(pdN_rand.lambda);
            end
            NP_E(j,i) = std(NP_rand);
            Broadening_E(j,i) = std(Broadening_rand);
            
            h = histogram(mydata,15,'Normalization','pdf',...
                'FaceColor',rgb('DodgerBlue'),'LineWidth',2,'FaceAlpha',0.7);
            xlabel(sprintf('counts in %.2f sec',mean(sstime)));
            ylabel('Frequency');
            PrettyFigureFormat
            % Gaussian PDF
            pdfG = @(b) 1/sqrt(2*pi*pdG.sigma.^2) .* exp(-(b-pdG.mu).^2/(2*pdG.sigma.^2));
            % Poisson PDF
            pdfP = @(b) 1/sqrt(2*pi*pdN.lambda) .* exp(-(b-pdN.lambda).^2/(2*pdN.lambda));
            % Convolution PDF
            %pdfC = @(b) Convolution(p(1),p(2),p(3),p(4),b);
            %pdfC = @(b) Convolution(Conv.lambda,Conv.sigma,Conv.mu,Conv.par4,b);
            hold on
            %b = (h.BinEdges(1:end-1)+h.BinWidth/2);
            b = linspace(min(mydata),max(mydata),100);
            g=plot(b,pdfG(b),'Color',rgb('IndianRed'),'LineWidth',4);
            p=plot(b,pdfP(b),'Color',rgb('GoldenRod'),'LineWidth',4);
            %c=plot(b,pdfC(b),'Color',rgb('Blue'),'LineWidth',4);
            leg=legend([h g p],...
                sprintf('%.0f scans - PSR %.0f',numel(R.RunList),i),...
                sprintf('Gaussian'),...
                sprintf('Poisson\nNP=%.1f\\pm%.1f\n\\sigma_{eq}=(%.2f\\pm%.2f)mV',pdG.sigma/sqrt(pdN.lambda),NP_E(j,i),Broadening(j,i),Broadening_E(j,i)),...
                ...%sprintf('Convolution'),...
                'location','northwest');
            leg.Color = 'none'; legend boxoff;
            hold off
            PrettyFigureFormat
            set(gca,'FontSize',12);
            disp(pdG.sigma/sqrt(pdN.lambda));

        end
        %export_fig(fig1,savefile1);

    end
    
    save('SamakKNM2_NPBroadeningsInRW123PSR1234_mV.mat','NP','Broadening','NP_E','Broadening_E');
    
elseif strcmp(Mode,'Periodwise')
        %% Rate Evolution --> mV equivalent
    myMainTitle = sprintf('KATRIN - KNM2 - FPD @E_0-300eV - Non-Poissonian Components');
    maintitle   = myMainTitle;
    fig1      = figure('Name',sprintf('KATRIN - KNM2 Scanwise Background'),...
        'NumberTitle','off','rend','painters','pos',[10 10 1400 1000]);
    a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
    a.FontSize=24;a.FontWeight='bold';

    savefile1    = sprintf('plots/KNM2_RM300_EffectivePotentialFluctuation.0f_PseudoRings_NP.png');


    % Loop on Period
    for j=1:3


        RunList   = ['KNM2_RW' num2str(j)];

        % Read Data
        DataType  = 'Real';
        RunAnaArg = {'RunList',RunList,'DataType',DataType,...
            'FSDFlag','BlindingKNM2','ELossFlag','KatrinT2',...
            'AnaFlag','StackPixel','RingMerge','Full','NonPoissonScaleFactor',1};
        R        = MultiRunAnalysis(RunAnaArg{:});
        range = 40;               % fit range in eV below endpoint        
        R.exclDataStart = R.GetexclDataStart(range); % find correct data, where to cut spectrum
        R.ROIFlag=ROIstr; R.SetROI;

        % Slow Control Data
        p1 = (R.SingleRunData.WGTS_MolFrac_TT'+0.5*R.SingleRunData.WGTS_MolFrac_HT'+0.5*R.SingleRunData.WGTS_MolFrac_DT')./mean((R.SingleRunData.WGTS_MolFrac_TT'+0.5*R.SingleRunData.WGTS_MolFrac_HT'+0.5*R.SingleRunData.WGTS_MolFrac_DT')).*R.SingleRunData.WGTS_CD_MolPerCm2'./mean(R.SingleRunData.WGTS_CD_MolPerCm2');
        p2 = mean(R.SingleRunData.qU_RM,1); p2=p2-mean(p2);

        %% HV Drift Correction
        FirstDayPeriod1 = datetime('02-Oct-2019 14:52:19');
        TimeLineDaysFirstDayPeriod1 = days(R.SingleRunData.StartTimeStamp-FirstDayPeriod1);
        if strcmp(HVDrift,'ON')
            HVdriftPerPixel = 1.5*(TimeLineDaysFirstDayPeriod1) * 6.3e-3;
        else
            HVdriftPerPixel = 0;
        end
        if strcmp(Corr,'Fabian')
            R.RMCorrection('QAplots',QAplots);
        end
        count  = R.SingleRunData.TBDIS_RM;
        sstime = mean(R.SingleRunData.qUfrac_RM,1).*R.SingleRunData.TimeSec;
        rate   = count./sstime;
        if strcmp(Detrend,'ON')
            rate = detrend(rate)+mean(rate);
        end
        if strcmp(Corr,'Fabian')
          corrcount_norm = rate .* mean(sstime);
        elseif strcmp(Corr,'Thierry')
          cf     = R.RMRateErosCorrectionqUActivity;
          corrcount_norm = rate.*cf .* mean(sstime);
        end

        subplot(1,3,j);

        %% Stacked Pixel: Non-Poissonian Component?
        lambda=2;
        sigma=1;
        mu=0.1;
        par4=1;
        x=-2:0.0001:2;
        p = [lambda,sigma,mu,par4,x];
        [N,edges] = histcounts(corrcount_norm,15,'Normalization','probability');
        means = zeros(1,numel(N));
        for n=1:(numel(edges)-1)
            means(n) = edges(n)+(edges(n+1)-edges(n))./2;
        end
        means = means - mean(edges);
        Poisson = @(lambda,x) lambda.^x.*exp(1).^(-lambda)./gamma(x+1);
        Gauss   = @(sigma,mu,x) exp(-0.5.*((x-mu)./sigma).^2)./(sigma.*sqrt(2.*pi));
        Convolution = @(lambda,sigma,mu,par4,x) par4.*conv(Poisson(lambda,x),Gauss(sigma,mu,x),'same');
        %Conv = fit(means',N',stairs(Convolution),'StartPoint',p,'lower',[1e5,0,1e5,0],'upper',[1e8,1e2,1e8,10]);
        %plot(Conv);

        pdG = fitdist(corrcount_norm','Normal');
        pdN = fitdist(corrcount_norm','poisson');
        mydata=corrcount_norm(abs(corrcount_norm-pdG.mu)<(3*pdG.sigma))';
        pdG = fitdist(mydata,'Normal');
        pdN = fitdist(mydata,'poisson'); 
        Broadening(j) = sqrt(pdG.sigma^2-pdN.lambda)/(mean(sstime)*737.8) * 1e3 * 117 / numel(R.PixList);
        NP(j)         = pdG.sigma/sqrt(pdN.lambda);
        
        %% Error calculation: Bootstrap
        mydata_rand=zeros(numel(mydata),1);
        Broadening_rand=zeros(NBoot,1);
        NP_rand=zeros(NBoot,1);
        for l=1:NBoot
              for k=1:numel(mydata)
                  mydata_rand(k) = mydata(randi(numel(mydata)));
              end
              pdG_rand = fitdist(mydata_rand,'Normal');
              pdN_rand = fitdist(mydata_rand,'poisson');
              Broadening_rand(l) = sqrt(abs(pdG_rand.sigma^2-pdN_rand.lambda))/(mean(sstime)*737.8) * 1e3 * 117 / numel(R.PixList);
              NP_rand(l)         = pdG_rand.sigma/sqrt(pdN_rand.lambda);
        end
        NP_E(j) = std(NP_rand);
        Broadening_E(j) = std(Broadening_rand);

        h = histogram(mydata,15,'Normalization','pdf',...
            'FaceColor',rgb('DodgerBlue'),'LineWidth',2,'FaceAlpha',0.7);
        xlabel(sprintf('counts in %.2f sec',mean(sstime)));
        ylabel('Frequency');
        PrettyFigureFormat
        % Gaussian PDF
        pdfG = @(b) 1/sqrt(2*pi*pdG.sigma.^2) .* exp(-(b-pdG.mu).^2/(2*pdG.sigma.^2));
        % Poisson PDF
        pdfP = @(b) 1/sqrt(2*pi*pdN.lambda) .* exp(-(b-pdN.lambda).^2/(2*pdN.lambda));
        % Convolution PDF
        %pdfC = @(b) Convolution(p(1),p(2),p(3),p(4),b);
        %pdfC = @(b) Convolution(Conv.lambda,Conv.sigma,Conv.mu,Conv.par4,b);
        hold on
        %b = (h.BinEdges(1:end-1)+h.BinWidth/2);
        b = linspace(min(mydata),max(mydata),100);
        g=plot(b,pdfG(b),'Color',rgb('IndianRed'),'LineWidth',4);
        p=plot(b,pdfP(b),'Color',rgb('GoldenRod'),'LineWidth',4);
        %c=plot(b,pdfC(b),'Color',rgb('Blue'),'LineWidth',4);
        leg=legend([h g p],...
            sprintf('RW%i',j),...
            sprintf('Gaussian'),...
            sprintf('Poisson\nNP=%.1f\\pm%.1f\n\\sigma_{eq}=(%.2f\\pm%.2f)mV',pdG.sigma/sqrt(pdN.lambda),NP_E(j),Broadening(j),Broadening_E(j)),...
            ...%sprintf('Convolution'),...
            'location','northwest');
        leg.Color = 'none'; legend boxoff;
        hold off
        PrettyFigureFormat
        set(gca,'FontSize',12);
        disp(pdG.sigma/sqrt(pdN.lambda));

    end
    
    save('SamakKNM2_NPBroadeningsInRW123uniform_mV.mat','NP','Broadening','NP_E','Broadening_E');
    
elseif strcmp(Mode,'Uniform')                                                                   %Needs Rework!
            %% Rate Evolution --> mV equivalent
    myMainTitle = sprintf('KATRIN - KNM2 - FPD @E_0-300eV - Non-Poissonian Components');
    maintitle   = myMainTitle;
    fig1      = figure('Name',sprintf('KATRIN - KNM2 Scanwise Background'),...
        'NumberTitle','off','rend','painters','pos',[10 10 1400 1000]);
    a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
    a.FontSize=24;a.FontWeight='bold';

    savefile1    = sprintf('plots/KNM2_RM300_EffectivePotentialFluctuation.0f_PseudoRings_NP.png');


    % Loop on Period
    for j=1:3


        RunList   = ['KNM2_RW' num2str(j)];

        % Read Data
        DataType  = 'Real';
        RunAnaArg = {'RunList',RunList,'DataType',DataType,...
            'FSDFlag','BlindingKNM2','ELossFlag','KatrinT2',...
            'AnaFlag','StackPixel','RingMerge','Full','NonPoissonScaleFactor',1};
        R        = MultiRunAnalysis(RunAnaArg{:});
        range = 40;               % fit range in eV below endpoint        
        R.exclDataStart = R.GetexclDataStart(range); % find correct data, where to cut spectrum
        R.ROIFlag=ROIstr; R.SetROI;

        % Slow Control Data
        p1 = (R.SingleRunData.WGTS_MolFrac_TT'+0.5*R.SingleRunData.WGTS_MolFrac_HT'+0.5*R.SingleRunData.WGTS_MolFrac_DT')./mean((R.SingleRunData.WGTS_MolFrac_TT'+0.5*R.SingleRunData.WGTS_MolFrac_HT'+0.5*R.SingleRunData.WGTS_MolFrac_DT')).*R.SingleRunData.WGTS_CD_MolPerCm2'./mean(R.SingleRunData.WGTS_CD_MolPerCm2');
        p2 = mean(R.SingleRunData.qU_RM,1); p2=p2-mean(p2);

        %% HV Drift Correction
        FirstDayPeriod1 = datetime('02-Oct-2019 14:52:19');
        TimeLineDaysFirstDayPeriod1 = days(R.SingleRunData.StartTimeStamp-FirstDayPeriod1);
        if strcmp(HVDrift,'ON')
            HVdriftPerPixel = 1.5*(TimeLineDaysFirstDayPeriod1) * 6.3e-3;
        else
            HVdriftPerPixel = 0;
        end
        if strcmp(Corr,'Fabian')
            R.RMCorrection('QAplots',QAplots);
        end
        count  = R.SingleRunData.TBDIS_RM;
        sstime = mean(R.SingleRunData.qUfrac_RM,1).*R.SingleRunData.TimeSec;
        rate   = count./sstime;
        if strcmp(Detrend,'ON')
            rate = detrend(rate)+mean(rate);
        end
        if strcmp(Corr,'Fabian')
          Corrcount_norm{j} = rate .* mean(sstime);
        elseif strcmp(Corr,'Thierry')
          cf     = R.RMRateErosCorrectionqUActivity;
          Corrcount_norm{j} = rate.*cf .* mean(sstime);
        end
    end

    %% Stacked Pixel: Non-Poissonian Component?
    corrcount_norm = [Corrcount_norm{1} Corrcount_norm{2} Corrcount_norm{3}];
    lambda=2;
    sigma=1;
    mu=0.1;
    par4=1;
    x=-2:0.0001:2;
    p = [lambda,sigma,mu,par4,x];
    [N,edges] = histcounts(corrcount_norm,15,'Normalization','probability');
    means = zeros(1,numel(N));
    for n=1:(numel(edges)-1)
        means(n) = edges(n)+(edges(n+1)-edges(n))./2;
    end
    means = means - mean(edges);
    Poisson = @(lambda,x) lambda.^x.*exp(1).^(-lambda)./gamma(x+1);
    Gauss   = @(sigma,mu,x) exp(-0.5.*((x-mu)./sigma).^2)./(sigma.*sqrt(2.*pi));
    Convolution = @(lambda,sigma,mu,par4,x) par4.*conv(Poisson(lambda,x),Gauss(sigma,mu,x),'same');
    %Conv = fit(means',N',stairs(Convolution),'StartPoint',p,'lower',[1e5,0,1e5,0],'upper',[1e8,1e2,1e8,10]);
    %plot(Conv);

    pdG = fitdist(corrcount_norm','Normal');
    pdN = fitdist(corrcount_norm','poisson');
    mydata=corrcount_norm(abs(corrcount_norm-pdG.mu)<(3*pdG.sigma))';
    pdG = fitdist(mydata,'Normal');
    pdN = fitdist(mydata,'poisson'); 
    Broadening = sqrt(pdG.sigma^2-pdN.lambda)/(mean(sstime)*737.8) * 1e3 * 117 / numel(R.PixList);

    h = histogram(mydata,15,'Normalization','pdf',...
        'FaceColor',rgb('DodgerBlue'),'LineWidth',2,'FaceAlpha',0.7);
    xlabel(sprintf('counts in %.2f sec',mean(sstime)));
    ylabel('Frequency');
    PrettyFigureFormat
    % Gaussian PDF
    pdfG = @(b) 1/sqrt(2*pi*pdG.sigma.^2) .* exp(-(b-pdG.mu).^2/(2*pdG.sigma.^2));
    % Poisson PDF
    pdfP = @(b) 1/sqrt(2*pi*pdN.lambda) .* exp(-(b-pdN.lambda).^2/(2*pdN.lambda));
    % Convolution PDF
    %pdfC = @(b) Convolution(p(1),p(2),p(3),p(4),b);
    %pdfC = @(b) Convolution(Conv.lambda,Conv.sigma,Conv.mu,Conv.par4,b);
    hold on
    %b = (h.BinEdges(1:end-1)+h.BinWidth/2);
    b = linspace(min(mydata),max(mydata),100);
    g=plot(b,pdfG(b),'Color',rgb('IndianRed'),'LineWidth',4);
    p=plot(b,pdfP(b),'Color',rgb('GoldenRod'),'LineWidth',4);
    %c=plot(b,pdfC(b),'Color',rgb('Blue'),'LineWidth',4);
    leg=legend([h g p],...
        sprintf('%.0f scans',numel(R.RunList)),...
        sprintf('Gaussian'),...
        sprintf('Poisson \n NP = %.1f \n \\sigma_{eq} = %.2f mV',pdG.sigma/sqrt(pdN.lambda),Broadening),...
        ...%sprintf('Convolution'),...
        'location','northwest');
    leg.Color = 'none'; legend boxoff;
    hold off
    PrettyFigureFormat
    set(gca,'FontSize',12);
    disp(pdG.sigma/sqrt(pdN.lambda));
    ylim([0 10e-5]);
    
    save('SamakKNM2_NPBroadeningsInRW123uniform_mV.mat','NP','Broadenings');
end