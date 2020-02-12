function [exE_broadened,exP_broadened]= FSD_Shift(exE,exP,sigma,varargin)
% ------------------------------------------------------------------------------------------------
% for sigma==0
% input:
% exE =  excitation energy
% exP =  excitation probability
% sigma =  width of gaussians (e.g. sigma doppler effect)
% output:
% excitation energy and excitation probability of broadened FSD
%
% new method: do not use binned gaussians but acutal functions!
% Lisa Schlüter (25th Nov, 2019)
% ------------------------------------------------------------------------------------------------

p=inputParser;
p.addParameter('SanityPlot','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('MultiGauss_RelPos','',@(x) isfloat(x));  %3 gaussians instead of using 1 gaussian per energy (for 3 RW settings)
p.addParameter('MultiGauss_Weights','',@(x) isfloat(x)); %2 gaussians instead of using 1 gaussian per energy 
p.parse(varargin{:});
SanityPlot         = p.Results.SanityPlot;
MultiGauss_RelPos  = p.Results.MultiGauss_RelPos;
MultiGauss_Weights = p.Results.MultiGauss_Weights;

% construct gaussians: every point is gaussian with width sigmaDE
gauss = @(x,mu,sigma) 1/(sigma*sqrt(2*pi))*exp(-0.5*(x-mu).^2./sigma.^2);
if isempty(MultiGauss_RelPos)
    exP_SingleGaussians = @(x) exP'.*gauss(x,exE',sigma);
else
    exP_SingleGaussians = @(x) MultiGauss_Weights(1).*exP'.*gauss(x,exE'+MultiGauss_RelPos(1),sigma)+ ...
                               MultiGauss_Weights(2).*exP'.*gauss(x,exE'+MultiGauss_RelPos(2),sigma)+...
                               MultiGauss_Weights(3).*exP'.*gauss(x,exE'+MultiGauss_RelPos(3),sigma);
end
exP_SumGaussians = @(x) sum(exP_SingleGaussians(x));

% rebinning: define integration boundaries
Elow = zeros(numel(exE)+1,1);  % lower boundary
Elow(1) = -inf;
Elow(2:end) = exE(1:end);
Eup = zeros(numel(exE)+1,1);   % upper boundary
Eup(1:end-1) = exE;
Eup(end) = inf;

% rebinned energies: take center of each bin
exE_rebin = Elow+(Eup-Elow)/2;
exE_rebin(1)   = exE(1)-(exE(2)-exE(1))/2;         % first bin (-inf to exE(1))
exE_rebin(end) = exE(end)+(exE(end)-exE(end-1))/2; % last bin ( exE(end) to inf)

% rebin: integrate gaussians
exP_broadened = arrayfun(@(a,b) integral(exP_SumGaussians,a,b),Elow,Eup);

% keep only entries with non-zero probability
keepIndex     = exP_broadened~=0;
exE_broadened = exE_rebin(keepIndex)';
exP_broadened = exP_broadened(keepIndex)';

% plot
if strcmp(SanityPlot,'ON')
    pFSD = plot(exE,exP,'-','LineWidth',2);
    hold on;
    pGaus_rebin = plot(exE_broadened,exP_broadened,'-','LineWidth',pFSD.LineWidth);
    hold off;
    legend([pFSD,pGaus_rebin],'Regular FSD','Broadened FSD (rebinned)');
    legend boxoff;
    PrettyFigureFormat('FontSize',24);
    xlim([-1 60])
    ylim([0,max(exP)*1.1])
end
end