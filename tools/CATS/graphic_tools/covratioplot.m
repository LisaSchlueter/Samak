function covratioplot(s)
% Change in covariance. COVRATIO plot of the fit performed by CATS
%
% covratioplot(s) plots the covratio of the fit performed by CATS where s is
% the output result of CATS ( s = cats(...) ).
%
% Copyright 2005-2011 Guillaume MENTION, CEA Saclay
% $Revision: 0.99$  $Date: 2011/11/15$

nobs = s.nobs;
nsys = s.nsys;
% npar = s.npar;

clf;
barh(s.covratio(1:nobs));
if nsys==0
    colormap('default')
else
    colormap(lines(nsys));
end
cmap = [.5*ones(1,3); circshift(lines(nsys),-1)];
hold on;
for i=1:nsys
    barh(nobs+i,s.covratio(nobs+i),'FaceColor',cmap(i,:));
end
hold off;
box on;
axis([0 max(s.covratio)*1.2 0 length(s.covratio)+1 ]);
grid off;
% set(gca,'Xtick',round(linspace(1,length(s.covratio),min(10,nobs+nsys))));
if isempty(s.names)
    set(gca,'YDir','Reverse','Ytick',1:(nobs+nsys),'YtickLabel',1:(nobs+nsys));
    grid on;
    set(gca,'YGrid','off');
    PrettyFigureFormat;
else
    set(gca,'YDir','Reverse','Ytick',1:(nobs+nsys),'YtickLabel',s.names(1:(nobs+nsys)));
    grid on;
    set(gca,'YGrid','off');
    PrettyFigureFormat;
%     rotateticklabel(gca,45);
end

xlabel('Covariance ratios (  \chi^2_{min} / ndof ratios )');

line([1 1],[1 length(s.covratio)],'LineWidth',2,'LineStyle','--',...
    'Color',.5*ones(3,1));
if (max(s.covratio>5))
    set(gca,'XLim',[0 5]);
    
end
