function [exE_broadened,exP_broadened]= FSD_Conv(exE,exP,sigma,varargin)
% ------------------------------------------------------------------------------------------------
% broaden final state distribution: replace every discrete state with a gaussian of width sigma
% input:
% exE =  excitation energy
% exP =  excitation probability
% sigma =  width of gaussians (e.g. sigma doppler effect)
% output: 
% excitation energy and excitation probability of broadened FSD
%
% Lisa Schlüter (Nov, 2019)
% ------------------------------------------------------------------------------------------------
p=inputParser;
p.addParameter('SanityPlot','OFF',@(x)ismember(x,{'ON','OFF'}));
p.parse(varargin{:});
SanityPlot = p.Results.SanityPlot;

exE = exE';
exP = exP';
% use fine energy binning for gaussians
exE_fb = (min(exE)-10):0.0001:(max(exE)+10); % energy vector: finer binning

% construct gaussians: every point is gaussian with width sigmaDE
exP_SingleGaussians = exP.*gaussian(exE_fb,exE,sigma)./simpsons(repmat(exE_fb,[numel(exE),1])',gaussian(exE_fb,exE,sigma)')'; %
exP_SumGaussians = sum(exP_SingleGaussians,1);

% rebin back to original energy vector
exP_broadened = zeros(size(exE));
nfb = zeros(numel(exE),1);
for i=1:numel(exE)-1
    e_tmp = exE_fb(exE_fb>= exE(i) & exE_fb <exE(i+1));%energies within the fine binning
    p_tmp = exP_SumGaussians(exE_fb>= exE(i) & exE_fb <exE(i+1)); %probabilities within the fine binning
    nfb(i) = numel(e_tmp);
    if all(p_tmp==0) %ignore entries with only zero -> they don't have an impact on TBDDS
    elseif numel(e_tmp)==1
        exP_broadened(i) =p_tmp.*(exE_fb(2)-exE_fb(1)); %if only 1 entry in bin -> no integration necessary
    else
        exP_broadened(i) = simpsons(e_tmp,p_tmp); %inteprate over bin
    end
end

% keep only entries with non-zero probability
keepIndex     = exP_broadened~=0;
exE_broadened = exE(keepIndex)';
exP_broadened = exP_broadened(keepIndex)';

% plot
if strcmp(SanityPlot,'ON')
    pFSD = plot(exE,exP,'-','LineWidth',2);
    hold on;
    pGaus = plot(exE_fb,exP_SumGaussians,'-','LineWidth',pFSD.LineWidth);
    hold on;
    pGaus_rebin = plot(exE_broadened,exP_broadened,'-','LineWidth',pFSD.LineWidth);
    hold off;
    legend([pFSD,pGaus,pGaus_rebin],'Regular FSD','Broadened FSD','Broadened FSD (rebinned)');
    legend boxoff;
    PrettyFigureFormat;
    xlim([-1 60])
    ylim([0,max(exP)*1.1])
end
end

