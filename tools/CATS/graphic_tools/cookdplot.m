function cookdplot(s)
% Cook's Distance plot of the fit performed by CATS
%
% cookdplot(s) plots the Cook's distance of the fit performed by CATS where s is
% the output result of CATS ( s = cats(...) ).
%
% Copyright 2005-2011 Guillaume MENTION, CEA Saclay
% $Revision: 0.99$  $Date: 2011/11/15$

nobs = s.nobs;
nsys = s.nsys;
% npar = s.npar;

clf;
barh(s.cookd(1:nobs),'facecolor',rgb('CadetBlue'));
if nsys==0
    colormap('default')
else
    colormap(lines(nsys));
end
cmap = [.5*ones(1,3); circshift(lines(nsys),-1)];
hold on;
for i=1:nsys
    barh(nobs+i,s.cookd(nobs+i),'FaceColor',cmap(i,:));
end
hold off;
box on;
axis([ 0 max(s.cookd)*1.2 0 length(s.cookd)+1]);
grid off;
% set(gca,'Xtick',round(linspace(1,length(s.cookd),min(10,nobs+nsys))));
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

xlabel('Cook''s distance');

line([1 1],[1 length(s.cookd)],'LineWidth',2,'LineStyle','--',...
    'Color',.5*ones(3,1));
line([4/nobs 4/nobs],[1 length(s.cookd)],'LineWidth',1,'LineStyle','-.',...
    'Color',.5*ones(3,1));
if(max(s.cookd)>5)
    set(gca,'XLim',[0 5]);
end

end




