%%------------------------------------------------------------------%%
% Script to get non-poisson factor from background distribution
% stacked FPD or ringwise
% Fit 1x Gauss and 1x Poisson distribution
% 'NonPoissonScaleFactor' = sigma(Gauss)/sigma(Poisson)
% Lisa, 16.10.2019
%%------------------------------------------------------------------%%
%% Settings
function [NPfactor, NPfactorErr]= GetNPfactor(obj,varargin)
p=inputParser;
p.addParameter('RecomputeFlag','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('SanityPlot','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('SlopeCorr','OFF',@(x)ismember(x,{'ON','OFF'}));  % only working for stacked FPD
p.parse(varargin{:});
RecomputeFlag = p.Results.RecomputeFlag;
SanityPlot    = p.Results.SanityPlot;
SlopeCorr     = p.Results.SlopeCorr;

if strcmp(obj.AnaFlag,'StackPixel')
    nPseudoRing = 1; % define number of pseudo rings
elseif strcmp(obj.AnaFlag,'Ring')
    nPseudoRing = numel(obj.RingList);
    SlopeCorr = 'OFF';
    fprintf('No background slope correction for rings - Set SlopeCorr OFF \n')
end

%% label
savedir = [getenv('SamakPath'),'tritium-data/NPfactor/'];
MakeDir(savedir);

if strcmp(obj.DataType,'Real')
    runname = obj.RunData.RunName;
else
    runname = extractBefore(obj.RunData.RunName,'_E0');
end

if strcmp(obj.AnaFlag,'StackPixel')
% <<<<<<< HEAD
%     savename = [savedir,sprintf('NPfactor_%s_UniformFPD_SlopeCorr%s.mat',obj.RunData.RunName,SlopeCorr)];
% elseif strcmp(obj.AnaFlag,'Ring')
%     savename = [savedir,sprintf('NPfactor_%s_Ring%s.mat',obj.RunData.RunName,obj.RingMerge)];
% =======
    savename = [savedir,sprintf('NPfactor_%s_UniformFPD_SlopeCorr%s.mat',runname,SlopeCorr)];
elseif strcmp(obj.AnaFlag,'Ring')
    savename = [savedir,sprintf('NPfactor_%s_Ring%s.mat',runname,obj.RingMerge)];
%>>>>>>> 2680d1afe0066666933d7782d14d252046b5598f
end
%% load of compute
if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename);
else
    fprintf('Calculate NP factor for %s \n',extractBefore(savename,'.mat'));
    %% init
    BkgIndex = find(obj.RunData.qU>obj.ModelObj.Q_i,1);
    nBkg = numel(obj.RunData.qU(BkgIndex:end,1));
    
    Bkg = zeros(numel(obj.RunList).*nBkg,nPseudoRing); %reshaped: nruns*nsubruns, nrings
    Bkg_Subrun = zeros(nBkg,numel(obj.RunList),nPseudoRing); % nsubruns, nruns, nrings
    dNorm = cell(nPseudoRing,1);
    dPois = cell(nPseudoRing,1);
    BkgSlope = cell(nPseudoRing,1);
    %% get background counts
    % and normalize to average time spent in background region
    
    if strcmp(obj.AnaFlag,'StackPixel')
        TimeperSubRun =   (obj.SingleRunData.qUfrac(BkgIndex:end,:).*obj.SingleRunData.TimeSec);
        MeanTimeperSubRun = mean(obj.SingleRunData.qUfrac(BkgIndex:end,:).*obj.SingleRunData.TimeSec,2);
        Bkg_Subrun =  obj.SingleRunData.TBDIS(BkgIndex:end,:).*(MeanTimeperSubRun./TimeperSubRun);
        Bkg = reshape(Bkg_Subrun,[nBkg*numel(obj.RunList),1]);
    else
        for i=1:nPseudoRing
            Bkg_Subrun(:,:,i) =  obj.SingleRunData.TBDIS(BkgIndex:end,:,i)./(obj.SingleRunData.qUfrac(BkgIndex:end,:,i)...
                .*obj.SingleRunData.TimeSec).*mean(obj.SingleRunData.qUfrac(BkgIndex:end,:,i).*obj.SingleRunData.TimeSec,2);
            Bkg(:,i) = reshape(Bkg_Subrun(:,:,i),[nBkg*numel(obj.RunList),1]);
        end
    end
    %% correct for increased background level
    
    % get background level from single run fits (uniform FPD)
    if strcmp(obj.AnaFlag,'StackPixel') && strcmp(SlopeCorr,'ON')
        knm2_FitBKGSlopeTime(obj);
        Bkg_fit = obj.SingleRun_FitResults.chi2Stat.B;
        Bkg_fitErr= obj.SingleRun_FitResults.chi2Stat.BErr;
        Bkg_fit_corr = zeros(numel(obj.RunList),nBkg); % corrected fit results -> should be flat!
        
        LiveTime   = seconds(obj.SingleRunData.StartTimeStamp-obj.SingleRunData.StartTimeStamp(1)); % time in hours with respect to first run
        MeanBkgfit = wmean(Bkg_fit,1./Bkg_fitErr);
        refTime = LiveTime(find(abs(Bkg_fit-MeanBkgfit)==min(abs(Bkg_fit-MeanBkgfit))));%mean(LiveTime);%LiveTime(end)/2;% ;
        LiveTime = LiveTime-refTime;
        
        % fit linear slope
        [coeff,coeffErr,chi2,dof]= linFit(LiveTime',Bkg_fit,Bkg_fitErr);
        BkgSlopeCpsPerS = coeff(1);
        fprintf('Background Slope = %.4g mcps/day \n',BkgSlopeCpsPerS*60*60*24*1e3);
        
        
        if strcmp(SanityPlot,'ON')
            f1 = figure('Units','normalized','Position',[0.1,0.1,0.8,0.5]);
            e1 =errorbar(LiveTime./(60*60),Bkg_fit*1e3,Bkg_fitErr*1e3,'o','LineWidth',1.5,'Color',obj.PlotColor,'MarkerFaceColor',rgb('SkyBlue'));
            e1.CapSize = 0;
            hold on;
            l1 = plot(LiveTime./(60*60),(coeff(1).*LiveTime./(60*60)+coeff(2))*1e3,'LineWidth',3,'Color',rgb('Orange'));
            PrettyFigureFormat('FontSize',24);
            xlabel('Live time (h)');
            ylabel('Background rate (mcps)')
            xlim([min(LiveTime./(60*60))-20,max(LiveTime./(60*60))+20]);
            ylim([175,270]);
            leg = legend([e1,l1],'Scanwise fits (statistics only)',sprintf('Background slope = (%.2g \\pm %.1g) mcps/day \n',...
                BkgSlopeCpsPerS*60*60*24*1e3,coeffErr(1)*60*60*24*1e3));
            legend boxoff
            leg.Location = 'northwest';
            savedirplot = [savedir,'plots/'];
            MakeDir(savedirplot);
            saveplot = [savedirplot,sprintf('BkgSlope_%s.pdf',obj.RunData.RunName)];
            export_fig(f1,saveplot,'-painters');
        end
        
        % correct background for slope
        if strcmp(obj.AnaFlag,'StackPixel')
            Bkg_Subrun_corr = Bkg_Subrun-(BkgSlopeCpsPerS.*LiveTime.*MeanTimeperSubRun);
            Bkg_corr = reshape(Bkg_Subrun_corr,[nBkg*numel(obj.RunList),1]);
            CorrFactor = sum(Bkg)/sum(Bkg_corr);
            Bkg_corr = Bkg_corr.*CorrFactor;
        elseif strcmp(obj.AnaFlag,'Ring')
            
            Bkg_Subrun_corr = zeros(nBkg,numel(obj.RunList),nPseudoRing);
            for i=1:nPseudoRing
                BkgSlopeCpsPerS_Ring = BkgSlopeCpsPerS/numel(B.RingPixList{i});
                Bkg_Subrun_corr(:,:,i) = Bkg_Subrun(:,:,i)-(BkgSlopeCpsPerS_Ring.*LiveTime.*MeanTimeperSubRun);
            end
            Bkg_corr = reshape(Bkg_Subrun_corr,[nBkg*numel(obj.RunList)*nPseudoRin,1]);
            CorrFactor = sum(Bkg)/sum(Bkg_corr);
            Bkg_corr = Bkg_corr.*CorrFactor;
            fprintf('slope correcton only for uniform FPD \n')
        end
        
        if strcmp(SanityPlot,'ON') && strcmp(SlopeCorr,'ON')
            % fit linear slope to uncorr
            [coeff_uncorr,coeffErr_uncorr,~,~]= linFit(LiveTime',mean(Bkg_Subrun)',sqrt(mean(Bkg_Subrun))');
            [coeff1,coeffErr1,~,~]= linFit(LiveTime',mean(Bkg_Subrun_corr)',sqrt(mean(Bkg_Subrun_corr))');
            
            f1 = figure('Units','normalized','Position',[0.1,0.1,0.8,0.5]);
            e_uncorr =errorbar(LiveTime./(60*60),mean(Bkg_Subrun),sqrt(mean(Bkg_Subrun)),'o','LineWidth',1.5,'Color',rgb('Silver'),...
                'MarkerFaceColor',rgb('Silver'));
            e_uncorr.CapSize = 0;
            hold on;
            l_uncorr = plot(LiveTime./(60*60),coeff_uncorr(1)*LiveTime+coeff_uncorr(2),'Color',rgb('Silver'),'LineWidth',1.5);
            e1 =errorbar(LiveTime./(60*60),mean(Bkg_Subrun_corr),sqrt(mean(Bkg_Subrun_corr)),'o','LineWidth',1.5,'Color',obj.PlotColor,'MarkerFaceColor',rgb('SkyBlue'));
            e1.CapSize = 0;
            l_corr = plot(LiveTime./(60*60),coeff1(1)*LiveTime+coeff1(2),'Color',e1.Color,'LineWidth',1.5);
            
            PrettyFigureFormat('FontSize',24);
            xlabel('Live time (h)');
            ylabel('Background per HV setpoint (counts)')
            xlim([min(LiveTime./(60*60))-20,max(LiveTime./(60*60))+20]);
            
            leg = legend([e_uncorr,e1],'Uncorrected', 'Slope correction');
            legend boxoff
            leg.Location = 'northwest';
            savedirplot = [savedir,'plots/'];
            MakeDir(savedirplot);
            saveplot = [savedirplot,sprintf('BkgSlope_SlopeCorrection_%s.pdf',obj.RunData.RunName)];
            export_fig(f1,saveplot);
        end
    end
    %% fit distributions
    for i=1:nPseudoRing
        switch SlopeCorr
            case 'ON'
                dPois{i} = fitdist(Bkg_corr(:,i),'Poisson');
                dNorm{i} = fitdist(Bkg_corr(:,i),'Normal');
            case 'OFF'
                dPois{i} = fitdist(Bkg(:,i),'Poisson');
                dNorm{i} = fitdist(Bkg(:,i),'Normal');
        end
    end
    %% fit parameter and uncertainties
    SigmaGauss    = cellfun(@(x) x.sigma,dNorm); % sigma of gaussian distirbution
    SigmaErrGauss = cellfun(@(x) sqrt(x.ParameterCovariance(2,2)),dNorm); % 1 sigma uncertainty on sigma(Gauss)
    
    SigmaPoiss    = cellfun(@(x) sqrt(x.lambda),dPois); % sigma of poisson distribution
    SigmaErrPoiss = cellfun(@(x) x.ParameterCovariance^(1/4),dPois); % sqrt(sqrt(dPois{ring}.ParameterCovariance)); % 1 sigma uncertainty on sigma(Poiss)
    
    NPfactor    = (SigmaGauss./SigmaPoiss)';
    
    %% bootstrap: 
    % random number from 1 to nSamples
    nSamples = 1e4;
    RandIndex = zeros(nPseudoRing,nSamples);
    RandIndex = randi([1 size(Bkg,1)],size(Bkg,1),nSamples);
    
    Bkg_Samples = zeros(nPseudoRing,size(Bkg,1),nSamples);
    NPfactor_Samples = zeros(nPseudoRing,nSamples);
    
    for i=1:nSamples
        for r=1:nPseudoRing
            progressbar(i/nSamples);
            Bkg_Samples(r,:,i) = Bkg(RandIndex(:,i),r);
            dPois_samples = fitdist(squeeze(Bkg_Samples(r,:,i))','Poisson');
            dNorm_samples = fitdist(squeeze(Bkg_Samples(r,:,i))','Normal');
            NPfactor_Samples(r,i) = dNorm_samples.sigma/sqrt(dPois_samples.lambda);          
        end
    end
    
    NPfactorErr = std(NPfactor_Samples,0,2);%sqrt((1./SigmaPoiss.*SigmaErrGauss).^2+(SigmaGauss./SigmaPoiss.^2.*SigmaErrPoiss).^2)';
    
    %% save
    save(savename,'NPfactor','NPfactorErr','SigmaGauss','SigmaErrGauss',...
        'SigmaPoiss','SigmaErrPoiss','Bkg','dPois','dNorm','NPfactor_Samples','Bkg_Samples');
    
    if strcmp(SlopeCorr,'ON')
       save(savename,'Bkg_corr','-append');
    end
end

%% plot
if strcmp(SanityPlot,'ON')
    if strcmp(obj.AnaFlag,'Ring')
        ring = 1;
        dataleg = sprintf('Data (pseudo ring %.0f)',ring);
    else
        ring = 1;
        dataleg = 'Data';
    end
    
    
    f2 = figure('Units','normalized','Position',[0.1,0.1,0.6,0.6]);
    switch SlopeCorr
        case 'OFF'
            h1 = histogram(Bkg(:,ring),'Normalization','pdf','FaceColor',rgb('SteelBlue'),'FaceAlpha',0.8);
        case 'ON'
            h1 = histogram(Bkg_corr(:,ring),'Normalization','pdf','FaceColor',rgb('SteelBlue'),'FaceAlpha',0.8);
    end
    hold on
    plotx    = linspace(h1.BinEdges(1)-2,h1.BinEdges(end)+2,100)';
    plotNorm = pdf(dNorm{ring},plotx);
    plotx2 = 0:0.1:round(max(plotx));
    plotPois = interp1(0:1:round(max(plotx)),poisspdf(0:1:round(max(plotx)),dPois{ring}.lambda),plotx2,'spline');
    lN = plot(plotx,plotNorm,'LineWidth',4,'Color',rgb('Crimson'));
    lP = plot(plotx2,plotPois,'LineWidth',lN.LineWidth,'Color',[0.929,0.694,0.125]);
    hold off
    PrettyFigureFormat('FontSize',24)
    xlabel('Background counts');
    ylabel('Probability density');
    leg = legend([h1,lN,lP],dataleg,sprintf('Gauss   \\sigma = (%.2f \\pm %.2f) counts',SigmaGauss(ring),SigmaErrGauss(ring)),...
        sprintf('Poisson \\sigma = (%.2f \\pm %.2f) counts',SigmaPoiss(ring),SigmaErrPoiss(ring)));
    leg.Title.String = sprintf('Non Poissonian factor = %.1f%% \\pm %.1f%%',[((NPfactor(ring)-1).*100)';NPfactorErr(ring)'.*100]);
    %leg.Title.String = sprintf('Non Poissonian factor = %.1f%% \\pm %.1f%% \n',[((NPfactor-1).*100)';NPfactorErr'.*100]);
    legend boxoff;
    xlim([h1.BinEdges(1),h1.BinEdges(end)]);
    title(sprintf('%s',obj.DataSet))
    savedirplot = [savedir,'plots/'];
    MakeDir(savedirplot);
       
    if strcmp(obj.AnaFlag,'StackPixel')
    saveplot = [savedirplot,sprintf('GetNPfactor_%s_SlopeCorr%s.pdf',obj.RunData.RunName,SlopeCorr)];
    else
       saveplot = [savedirplot,sprintf('GetNPfactor_%s_Ring%.0f_SlopeCorrOFF.pdf',obj.RunData.RunName,ring)];   
    end
    export_fig(f2,saveplot,'-painters');
    fprintf('Save plot to %s \n',saveplot);
    
    %% plot 2: plot ringwise NP-factor and fit a linear slope
    if strcmp(obj.AnaFlag,'Ring')
        x = 1:numel(obj.RingList);
        f3 = figure('Units','normalized','Position',[0.1,0.1,0.6,0.6]);
        
        % plot weighted mean
        MeanNP = wmean(NPfactor',1./NPfactorErr.^2);
        chi2NP = sum((NPfactor'-MeanNP).^2./NPfactorErr.^2);
        pval   = 1-chi2cdf(chi2NP,numel(NPfactor)-1); %consistency with weighted mean
        lm  = plot(linspace(min(x)-0.3,max(x)+0.3,10),MeanNP.*ones(10,1),'--','LineWidth',3,'Color',rgb('DarkGray'));
        
        hold on;
        
        % plot NP factors with error bars
        e1 = errorbar(x,NPfactor,NPfactorErr,'o','LineWidth',3,'Color',obj.PlotColor,'MarkerSize',0.1);
        e1.CapSize = 0;
        plot(x,NPfactor,'o','LineWidth',0.1,'Color',obj.PlotColor,'MarkerSize',8,'MarkerFaceColor',obj.PlotColor);
        
        % plot linear fit
        [linFitPar,linFitErr,chi2,dof]= linFit(x',NPfactor',NPfactorErr); % fit linear slope
        l1 = plot(x,linFitPar(1).*x+linFitPar(2),'-','Color',obj.PlotColor,'LineWidth',e1.LineWidth);
        
        % appearance, labels,...
        xlabel('Ring');
        xticks(obj.RingList); set(gca,'XMinorTick','off');
        if strcmp(obj.RingMerge,'Full')
            xticklabels({'1,2,3','4,5,6','7,8,9','10,11,12'})
        elseif strcmp(obj.RingMerge,'Half')
            xticklabels({'1,2,3,4,5,6','7,8,9,10,11,12'})
        end
        ylim([0.995*min(NPfactor-NPfactorErr'),1.005*max(NPfactor+NPfactorErr')]);
        xlim([min(x)-0.3, max(x)+0.3]);
        PrettyFigureFormat('FontSize',24);
        ylabel('Non Poisson factor');
        
        leg = legend([l1,lm],...
            sprintf(' Linear fit slope = %.1g \\pm %.1g per pseudo ring , p-value = %.2f',linFitPar(1),linFitErr(1),1-chi2cdf(chi2,dof)),...
            sprintf(' Weighted mean = %.2f , p-value = %.2f',MeanNP,pval));
        leg.EdgeColor = rgb('Silver');
        leg.Location='northeast';
        leg.FontSize = get(gca,'FontSize');
        saveplot = [savedirplot,sprintf('GetNPfactor_%s_Ring%.0f_RingWiseSlope.pdf',obj.RunData.RunName,ring)];  
        export_fig(f2,saveplot,'-painters');
        fprintf('Save plot to %s \n',saveplot);
    end
end
end

