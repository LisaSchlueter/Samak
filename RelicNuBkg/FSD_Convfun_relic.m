function [exE_broadened,exP_broadened]= FSD_Convfun_relic(exE,exP,sigma,varargin)
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
p.addParameter('mnu',0,@(x) isfloat(x));
p.addParameter('BinVec',0,@(x) isfloat(x) || isempty(x)); % reference bin vector
p.addParameter('RebinMode','Integral',@(x)ismember(x,{'Fast','Integral'})); % rebinning method. warning: 'fast' does not work welll for small sigma
p.addParameter('RecomputeFlag','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('filename','',@(x)ischar(x));

p.parse(varargin{:});

SanityPlot         = p.Results.SanityPlot;
ZoomPlot           = p.Results.ZoomPlot;
mnu                = p.Results.mnu;
BinVec             = p.Results.BinVec;
RebinMode          = p.Results.RebinMode;
RecomputeFlag      = p.Results.RecomputeFlag;
filename           = p.Results.filename;

Q = 18575;
MH3 = 2.809431917483801e9;  % Tritium mass eV, PRC77,055502(2008)
MHe3 = 2.809413380751987e9; % Helium-3 mass eV, PRC77,055502(2008)
me = 510.9989461e3;          % Electron Mass eV  PDG 2019

NBins = numel(BinVec);
savedir = [getenv('SamakPath'),'inputs/FSD/FSDconv/'];
savefile = strrep(strrep(extractAfter(filename,'inputs/FSD/'),'.txt',''),'.dat','');

savename = sprintf('%s%s_Relics_Sigma%.0fmeV_Binning%.0f.mat',...
savedir,savefile,sigma*1e3,NBins);

if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename);
else

    sigma(sigma<1e-2) = 1e-2; % otherwise num. integration problems 

    % replacement distribution
    %gaussian: replace every final state with a gaussian with a width sigma
    distfun = @(x,mu,sigma) 1/(sigma*sqrt(2*pi))*exp(-0.5*(x-mu).^2./sigma.^2);
    
    exE = mnu - exE;% + Q.*me./(me+MH3+MHe3);
    exP_SinglePeak = @(x) exP'.*distfun(x,exE',sigma);
    exP_SumPeak = @(x) sum(exP_SinglePeak(x));

    switch RebinMode
        case 'Integral' %proper integration
            %% rebinning step 1: define integration boundaries: keep the similar binning
            % lower boundary
            BinWidth = abs((BinVec(1)-BinVec(end))./(NBins-1));
            Elow = BinVec-(BinWidth./2); %(exE(end):BinWidth:(exE(1)-BinWidth));
            % upper boundary
            Eup = zeros(NBins,1);
            Eup(1:end-1) = BinVec(1:end-1)+(BinWidth./2); %(exE(end)+BinWidth:BinWidth:exE(1)); 
            Eup(end) = inf;

            %% rebinning step 2 probabilities: integrate gaussians
            %WayPoints = min(Elow):min(sigma):max(Emax);
            WayPoints = min(Elow(~isinf(Elow))):min(sigma):max(Eup(~isinf(Eup)));

            exP_rebin = arrayfun(@(a,b) integral(exP_SumPeak,a,b,'AbsTol',1e-9),Elow,Eup);
            exE_rebin = BinVec;

        case 'Fast' % fast and approximate
            exE_BinSize =  [exE(2)-exE(1),diff(exE)];
            exP_rebin = exP_SumPeak(exE).*exE_BinSize;
            exE_rebin = exE;

    %             % keep only entries with non-zero probability
    %             RmIndex     = exP_rebin~=0;
    %             exP_broadened{r} = exP_rebin(RmIndex);
    %             exE_broadened{r} = exE(RmIndex);
    end


    % keep only entries with non-zero probability (common number of bins)
    %exP_rebin = cell2mat(exP_rebin');
    %exE_rebin = cell2mat(exE_rebin');
    exP_broadened = exP_rebin';
    exE_broadened = exE_rebin';
    exE_broadened(isnan(exE_broadened(:,end)),end) = mean(exE_broadened(~isnan(exE_broadened(:,end)),end));

    MakeDir(savedir);
    save(savename,'exE_broadened','exP_broadened','exE','exP_SinglePeak','NBins');
end
%% plot
if strcmp(SanityPlot,'ON')
    exE_broadened_plot = exE_broadened;
    exP_broadened_plot = exP_broadened; 
    f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
    exE_plot = linspace(min(exE),max(exE),5000);
    pSingleGaus = plot(exE_plot,exP_SinglePeak(exE_plot),'-.','LineWidth',1.5,'Color',rgb('Silver'));
    hold on;
    pFSD = plot(exE,exP,'-','LineWidth',2,'Color',rgb('DodgerBlue'));
    pGaus_rebin = plot(exE_broadened_plot,exP_broadened_plot,'-','LineWidth',pFSD.LineWidth,'Color',rgb('Orange'));
    hold off;
    leg = legend([pFSD,pSingleGaus(1),pGaus_rebin],'Regular FSD','Replacement with Gaussians','Broadened FSD (rebinned)');
    %legend boxoff;
    leg.EdgeColor = rgb('LightGray');
    PrettyFigureFormat('FontSize',24);
    xlim([-1 40])
    xlabel('Excitation energy (eV)')
    ylabel('Probability')
    ylim([0,max(exP)*1.1])
    savedir = [getenv('SamakPath'),'inputs/Plasma/plots/'];
    MakeDir(savedir);
 
    savename = sprintf('%sFSD_%.3gSigma_Single.pdf',savedir,sigma);
    export_fig(f1,savename);
    
    fprintf('Save plot to %s\n',savename);
    
    if strcmp(ZoomPlot,'ON')
        pFSD.Marker = 'o';
        pFSD.MarkerSize = 8;
        pFSD.MarkerFaceColor = pFSD.Color;
        pGaus_rebin.Marker = 'o';
        pGaus_rebin.MarkerSize = 8;
        pGaus_rebin.MarkerFaceColor = pGaus_rebin.Color;
        xlim([23.1 24.1]);
        ylim([0 0.008]);
        
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