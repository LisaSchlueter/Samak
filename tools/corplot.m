function [ h cor c ] = corplot( cov , flag )
CC=squeeze(cov);
SS=sqrt(abs(diag(CC)));
cor=CC./(SS*SS');
n = size(cov,1)+1;
if nargin==2
    if flag
        hh = imagesc(cor);
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

% hold on
% [x,y] = meshgrid(1:numel(cov(1,:)),1:numel(cov(1,:)));
% tmp=round(reshape(cor,1,11*11),2);
% text(x(:)+0.5,y(:)+0.5,num2str(tmp'),'HorizontalAlignment','center','VerticalAlignment','middle');
% hold off

set(gca,'XAxisLocation','top');
set(gca,'XTick',.5+(1:(n-1)),'XTickLabel',1:(n-1));
set(gca,'YTick',.5+(1:(n-1)),'YTickLabel',1:(n-1));
set(gca,'TickLength',[0 0]);
caxis([-1 1]); cm=0:1/50:1; my_cm_blue=[1-cm;1-cm;ones(size(cm))]';my_cm_red=[ones(size(cm));1-cm;1-cm]';
colormap([flipud(my_cm_blue);[1 1 1];my_cm_red]);  
c = colorbar;

if nargout ~= 0
    h = hh;
end

end

