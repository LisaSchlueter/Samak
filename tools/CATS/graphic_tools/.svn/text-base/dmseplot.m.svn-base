function h = dmseplot(s)
% Change in Mean Squared Errors (chi2min/ndof).
%
% dmseplot(s) plots the mse change in the fit performed by CATS where s is
% the output result of CATS ( s = cats(...) ).
%
% Copyright 2005-2011 Guillaume MENTION, CEA Saclay
% $Revision: 0.99$  $Date: 2011/11/15$

nobs = s.nobs;
nsys = s.nsys;
% npar = s.npar;

clf;
h1 = barh(s.s2_i(1:nobs));
set(h1,'BaseValue',s.mse);
colormap(lines(nsys));
cmap = [.5*ones(1,3); circshift(lines(nsys),-1)];
hold on;
for i=1:nsys
    barh(nobs+i,s.s2_i(nobs+i),'FaceColor',cmap(i,:));
end
hold off;
box on;
axis([min(1,min(s.s2_i)/1.1) max(s.s2_i)*1.1 0 length(s.s2_i)+1 ]);
grid off;
% set(gca,'Xtick',round(linspace(1,length(s.s2_i),min(10,nobs+nsys))));
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

xlabel('change in mse (  \chi^2_{min} / ndof )');

line([1 1]*s.mse,[0 length(s.s2_i)+1],'LineWidth',2,'LineStyle','-',...
    'Color',.5*ones(3,1));
line([1 1]*s.mse*s.ndof/(s.ndof-1),...
    [0 length(s.s2_i)+1],'LineWidth',2,'LineStyle','--',...
    'Color',.5*ones(3,1));
line([1 1]*(s.mse*s.ndof-1)/(s.ndof-1),...
    [0 length(s.s2_i)+1],'LineWidth',2,'LineStyle',':',...
    'Color',.5*ones(3,1));

if (max(s.s2_i>5))
    set(gca,'XLim',[0 5]);
    
end
