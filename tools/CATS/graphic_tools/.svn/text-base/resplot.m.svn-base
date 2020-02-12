function resplot(s)
% Standardized residual plot of the fit performed by CATS
%
% resplot(s) plots the standardized residuals of the fit performed by CATS where s is
% the output result of CATS ( s = cats(...) ).
%
% Copyright 2005-2011 Guillaume MENTION, CEA Saclay
% $Revision: 0.99$  $Date: 2011/11/15$

nobs = s.nobs;
nsys = s.nsys;
% npar = s.npar;

clf;
barh(s.standres(1:nobs));
colormap(lines(nsys));
cmap = [.5*ones(1,3); circshift(lines(nsys),-1)];
hold on;
for i=1:nsys
    barh(nobs+i,s.standres(nobs+i),'FaceColor',cmap(i,:));
end
hold off;
box on;
axis([min(0,min(s.standres)*1.2) max(s.standres)*1.2 0 length(s.standres)+1 ]);
grid off;
% set(gca,'Xtick',round(linspace(1,length(s.leverage),min(10,nobs+nsys))));
if isempty(s.names)
    set(gca,'YDir','Reverse','Ytick',1:(nobs+nsys),'YtickLabel',1:(nobs+nsys));
    grid on;
    set(gca,'YGrid','off');
    PrettyFigureFormat;
else
    set(gca,'YDir','Reverse','Ytick',1:(nobs+nsys),'YtickLabel',s.names(1:(nobs+nsys)));
    PrettyFigureFormat;
    grid on;
    set(gca,'YGrid','off');
    %     rotateticklabel(gca,90);
    %     set(gca,'Visible','off');
end

xlabel('Standardized Residuals');

end