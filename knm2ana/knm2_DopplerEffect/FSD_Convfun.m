function [exE_broadened,exP_broadened]= FSD_Convfun(exE,exP,sigma,varargin)
% ------------------------------------------------------------------------------------------------
% broaden final state distribution: replace every discrete state with a gaussian of width sigma
% input:
% exE =  excitation energy
% exP =  excitation probability
% sigma =  width of gaussians (e.g. sigma doppler effect)
% output:
% excitation energy and excitation probability of broadened FSD (binned)
% excitation energies are rebinned using the weighted mean
%
% new method: do not use binned gaussians but analytical function
% Lisa Schlüter (25th Nov, 2019)
% ------------------------------------------------------------------------------------------------

p=inputParser;
p.addParameter('SanityPlot','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('ZoomPlot','OFF',@(x)ismember(x,{'ON','OFF'}));   % zoom into plot and save
p.addParameter('MultiPos','',@(x) isfloat(x) || isempty(x));     %3 gaussians instead of using 1 gaussian per energy (for 3 RW settings)
p.addParameter('MultiWeights','',@(x) isfloat(x) || isempty(x)); %3 gaussians instead of using 1 gaussian per energy
p.addParameter('Dist','Gauss',@(x)ismember(x,{'Gauss','Rect'}));
p.addParameter('BinningFactor',1,@(x) isfloat(x) || isempty(x)); % multiplies number of bins in rebinning
p.addParameter('RebinMode','Integral',@(x)ismember(x,{'Fast','Integral'})); % rebinning method. warning: 'fast' does not work welll for small sigma

p.parse(varargin{:});

SanityPlot         = p.Results.SanityPlot;
ZoomPlot           = p.Results.ZoomPlot;
MultiPos           = p.Results.MultiPos;
MultiWeights       = p.Results.MultiWeights;
BinningFactor      = p.Results.BinningFactor;
Dist               = p.Results.Dist;
RebinMode          = p.Results.RebinMode;

% replacement distribution
if strcmp(Dist,'Gauss')
    %gaussian: replace every final state with a gaussian with a width sigma
    distfun = @(x,mu,sigma) 1/(sigma*sqrt(2*pi))*exp(-0.5*(x-mu).^2./sigma.^2);
elseif strcmp(Dist,'Rect')
    % rectangle function: replace every final state with a linear distribution with a
    distfun = @(x,mu,sigma)  0+...%((x <  mu-sigma/2)  & (x >  mu+sigma/2)).*0 + ...   % zero for this condition
        ((x >= mu-sigma/2)  & (x <= mu+sigma/2)).*1/sigma;
end

nRings = ceil(numel(sigma)/3);% number of pseudo rings
exP_rebin = cell(nRings,1);
exE_rebin = cell(nRings,1);

for r=1:nRings
    if isempty(MultiPos)
        exP_SinglePeak = @(x) exP'.*distfun(x,exE',sigma(:,r));
    elseif nRings==1
        exP_SinglePeak = @(x) ( MultiWeights(1).*exP'.*distfun(x,exE'-MultiPos(1),sigma(1))+...
            MultiWeights(2).*exP'.*distfun(x,exE'-MultiPos(2),sigma(2))+...
            MultiWeights(3).*exP'.*distfun(x,exE'-MultiPos(3),sigma(3)));
    elseif size(MultiPos,1)==3
        exP_SinglePeak = @(x) ( MultiWeights(1).*exP'.*distfun(x,exE'-MultiPos(1,r),sigma(1,r))+...
            MultiWeights(2).*exP'.*distfun(x,exE'-MultiPos(2,r),sigma(2,r))+...
            MultiWeights(3).*exP'.*distfun(x,exE'-MultiPos(3,r),sigma(3,r)));
    else
        fprintf('MultiGauss only available for 3 gaussians \n');
    end
    exP_SumPeak = @(x) sum(exP_SinglePeak(x));
    
    switch RebinMode
        case 'Integral' %proper integration
            %% rebinning step 1: define integration boundaries: keep the similar binning
            % lower boundary
            Elow = zeros(numel(exE),1);
            Elow(1) = -inf;
            Elow(2:end) = mean([exE(1:end-1);exE(2:end)]);
            % upper boundary
            Eup = zeros(numel(exE),1);
            Eup(1:end-1) = mean([exE(1:end-1);exE(2:end)]);
            Eup(end) = inf;
            
            if ~isempty(BinningFactor)
                % enhance binning
                if BinningFactor>=2
                    RebinStop = floor(BinningFactor/2);
                    for i=1:RebinStop
                        Enew = Elow(2:end-1)+(Elow(2:end-1)-Eup(2:end-1))/2;
                        Elow = sort([Elow;Enew]);
                        Eup = sort([Eup;Enew]);
                    end
                end
            end
            %% rebinning step 2.1 probabilities: integrate gaussians
            exP_rebin{r} = arrayfun(@(a,b) integral(exP_SumPeak,a,b),Elow,Eup);
            
            %% rebinning step 2.2 energies: weighted mean of integration boundaries
            Weightfun = @(x) exP_SumPeak(x).*x;                                                       % weighted integral
            exE_rebin{r} = arrayfun(@(a,b) integral(Weightfun,a,b)./integral(exP_SumPeak,a,b),Elow,Eup); % weighted average
            %exE_rebin = mean([Elow';Eup'])'; % unweighted average -> replaced by weighted average
            %exE_rebin(end) = exE_rebin(end-1)+10;
                                                                   % exclude entries with zero probability
        case 'Fast' % fast and approximate
            exE_BinSize =  [exE(2)-exE(1),diff(exE)];
            exP_rebin{r} = exP_SumPeak(exE).*exE_BinSize;
            exE_rebin{r} = exE;
            
%             % keep only entries with non-zero probability
%             RmIndex     = exP_rebin~=0;
%             exP_broadened{r} = exP_rebin(RmIndex);
%             exE_broadened{r} = exE(RmIndex);
    end
    
end

% keep only entries with non-zero probability (common number of bins)
exP_rebin = cell2mat(exP_rebin');
exE_rebin = cell2mat(exE_rebin');
RmIndex   = exP_rebin==0;  
keepIndex = ~prod(RmIndex,2);
exP_broadened = exP_rebin(keepIndex,:)';
exE_broadened = exE_rebin(keepIndex,:)';
exE_broadened(isnan(exE_broadened(:,end)),end) = mean(exE_broadened(~isnan(exE_broadened(:,end)),end));
%% plot
if strcmp(SanityPlot,'ON')
    if nRings>1
        plotRing = 1; % plot for pseudo-ring number 1
        exE_broadened_plot = exE_broadened(plotRing,:);
        exP_broadened_plot = exP_broadened(plotRing,:);
    end
    f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
    exE_plot = linspace(min(exE),max(exE),5000);
    pSingleGaus = plot(exE_plot,exP_SinglePeak(exE_plot),'-.','LineWidth',1.5,'Color',rgb('Silver'));
    hold on;
    pFSD = plot(exE,exP,'-','LineWidth',2,'Color',rgb('DodgerBlue'));
    pGaus_rebin = plot(exE_broadened_plot,exP_broadened_plot,'-','LineWidth',pFSD.LineWidth,'Color',rgb('Orange'));
    hold off;
    switch Dist
        case 'Gauss'
            leg = legend([pFSD,pSingleGaus(1),pGaus_rebin],'Regular FSD','Replacement with Gaussians','Broadened FSD (rebinned)');
        case 'Rect'
            leg = legend([pFSD,pSingleGaus(1),pGaus_rebin],'Regular FSD','Replacement with rectangles','Broadened FSD (rebinned)');
    end
    %legend boxoff;
    leg.EdgeColor = rgb('LightGray');
    PrettyFigureFormat('FontSize',24);
    xlim([-1 60])
    xlabel('Excitation energy (eV)')
    ylabel('Probability')
    ylim([0,max(exP)*1.1])
    savedir = [getenv('SamakPath'),'inputs/Plasma/plots/'];
    MakeDir(savedir);
 
    if ~isempty(MultiPos)
        savename = sprintf('%sFSD_%.3gSigma_Multi%s.pdf',savedir,sigma,Dist);
        export_fig(f1,savename);
    else
        savename = sprintf('%sFSD_%.3gSigma_Single%s.pdf',savedir,igma,Dist);
        export_fig(f1,savename);
    end
    
    fprintf('Save plot to %s\n',savename);
    
    if strcmp(ZoomPlot,'ON')
        pFSD.Marker = 'o';
        pFSD.MarkerSize = 8;
        pFSD.MarkerFaceColor = pFSD.Color;
        pGaus_rebin.Marker = 'o';
        pGaus_rebin.MarkerSize = 8;
        pGaus_rebin.MarkerFaceColor = pGaus_rebin.Color;
        xlim([1 2.9])
        
        savenameZoom = strrep(savename,'.pdf','_Zoom.pdf');
        fprintf('Save zoomed plot to %s\n',savenameZoom);
        export_fig(f1,savenameZoom);
    end
end
end

%% old 
% OLD: normal mean
% exE_rebin = mean([Elow,Eup],2);
% exE_rebin(1)   = exE(1)-(exE(2)-exE(1))/2;         % first bin (-inf to exE(1))
% exE_rebin(end) = exE(end)+(exE(end)-exE(end-1))/2; % last bin ( exE(end) to inf)
% exE_broadened = exE_rebin(keepIndex)';  %exclude entries with zero probability