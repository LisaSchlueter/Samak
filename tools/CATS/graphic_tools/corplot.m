function [ h, cor ] = corplot( stats , flag )
% Correlation matrix plot from the covariance matrix
%
% corplot( cov , flag ) plots the correlation matrix. If flag is set to
% true, imagesc command is used for the plot, otherwise, pcolor is used
% instead.
%
% Copyright 2005-2010 Guillaume MENTION, CEA Saclay
% $Revision: 0.99$  $Date: 2011/01/28$
if isstruct(stats)
    cov = stats.covx;
else
    cov = stats;
end
CC=squeeze(cov);
SS=sqrt(diag(CC));
cor=CC./(SS*SS');
n = size(cov,1)+1;
if nargin==2
    if flag
        hh = imagesc(1.5:(n-.5),1.5:(n-.5),cor);
    else
        hh = pcolor(1:n,1:n,[cor zeros(n-1,1); zeros(1,n)]);
        axis ij;
        axis square;
    end
else
        hh = pcolor(1:n,1:n,[cor zeros(n-1,1); zeros(1,n)]);
        axis ij;
        axis square;
end
set(gca,'XAxisLocation','top');
if isstruct(stats)
    if isempty(stats.names)
        set(gca,'XTick',.5+(1:(n-1)),'XTickLabel',1:(n-1));
        set(gca,'YTick',.5+(1:(n-1)),'YTickLabel',1:(n-1));
        set(gca,'TickLength',[0 0]);
    else
        set(gca,'XTick',.5+(1:(n-1)),'XTickLabel',stats.names(stats.nobs+(1:(stats.nsys+stats.npar))));
        set(gca,'YTick',.5+(1:(n-1)),'YTickLabel',stats.names(stats.nobs+(1:(stats.nsys+stats.npar))));
        set(gca,'TickLength',[0 0]);
    end
else
    if n<=50
        set(gca,'XTick',.5+(1:(n-1)),'XTickLabel',1:(n-1));
        set(gca,'YTick',.5+(1:(n-1)),'YTickLabel',1:(n-1));
        set(gca,'TickLength',[0 0]);
    else
        set(gca,'XTick',[],'YTick',[]);
end

caxis([-1 1]);
% cm=0:1/30:1; my_cm_blue=[1-cm;1-cm;ones(size(cm))]';my_cm_red=[ones(size(cm));1-cm;1-cm]';
% colormap([flipud(my_cm_blue);[1 1 1];my_cm_red]);
colormap(flipud(lbmap(63,'RedBlue')));
colorbar;
% PrettyFigureFormat;
grid off;
if nargout ~= 0
    h = hh;
end

end

