function dfxsplot(s,flag)
% dfxs plot of the fit performed by CATS
%
% dfxsplot(s) plots the difference in fitted parameters of the fit performed by CATS where s is
% the output result of CATS ( s = cats(...) ).
%
% Copyright 2005-2011 Guillaume MENTION, CEA Saclay
% $Revision: 0.99$  $Date: 2011/11/15$

nobs = s.nobs;
nsys = s.nsys;
npar = s.npar;

clf;
if nargin==2
    if flag
        hh = imagesc(1:(nobs+nsys),1:(nsys+npar),s.dfxs);
    else
        hh = pcolor(1:(nobs+nsys),1:(nsys+npar),...
            [s.dfxs zeros(nsys+npar,1); zeros(1,nobs+nsys+1)]);
        axis ij;
        axis square;
    end
else
        hh = pcolor(1:(nobs+nsys+1),1:(nsys+npar+1),...
            [s.dfxs zeros(nsys+npar,1); zeros(1,nobs+nsys+1)]);
        axis ij;
        axis square;
end
set(gca,'XAxisLocation','top');
if isempty(s.names)
    set(gca,'XTick',.5+(1:(nobs+nsys)),'XTickLabel',1:(nobs+nsys));
    set(gca,'YTick',.5+(1:(nsys+npar)),'YTickLabel',1:(nsys+npar));
    set(gca,'TickLength',[0 0]);
else
    set(gca,'XTick',.5+(1:(nobs+nsys)),'XTickLabel',1:(nobs+nsys));
    set(gca,'YTick',.5+(1:(nsys+npar)),'YTickLabel',s.names(nobs+(1:(nsys+npar))));
    a=get(gca,'XTickLabel');
    %erase current tick labels from figure
    set(gca,'XTickLabel',[]);
    %get tick label positions
    b=get(gca,'XTick');
    c=get(gca,'YTick');
%     PrettyFigureFormat;
    th=text(b,repmat(c(1)-abs(c(2)-c(1)),nobs+nsys,1),s.names(1:(nobs+nsys)),'HorizontalAlignment','left',...
        'rotation',90,...
        'FontName',get(gca,'FontName'),...
        'FontSize',get(gca,'FontSize'),...
        'FontWeight',get(gca,'FontWeight'));
    set(gca,'TickLength',[0 0]);
end
mm=max(abs(s.dfxs(:)));
caxis([-mm mm]); cm=0:1/30:1; my_cm_blue=[1-cm;1-cm;ones(size(cm))]';my_cm_red=[ones(size(cm));1-cm;1-cm]';
colormap([flipud(my_cm_blue);[1 1 1];my_cm_red]);
colorbar;
set(gca,'DataAspectRatio',[1 1 1])
if ~isempty(s.names)
    text(mean(b),c(end)+4*abs(c(2)-c(1)),...
        'Standardized changes in fitted values\newline leaving out corresponding prior info',...
        'HorizontalAlignment','center',...
        'FontName',get(gca,'FontName'),...
        'FontSize',get(gca,'FontSize'),...
        'FontWeight',get(gca,'FontWeight'));
end
if nargout ~= 0
    h = hh;
end

end

