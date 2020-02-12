function h = dfxiplot(s,i)
% dfxiplot is a bar graph of change in fitted value of a given variable of the fit
% performed by CATS
%
% dfxiplot(s,i) uses the fit performed by CATS where s is
% the output result of CATS ( s = cats(...) ). By default if i is not
% specified this routine plots the standardized change in the last fitted
% parameter. Otherwise i could be an integer or a valid variable name to
% identify a given fitter parameter in the CATS structure. This routine
% plots a bar graph of a row of s.dfxs output.
%
% Copyright 2005-2011 Guillaume MENTION, CEA Saclay
% $Revision: 0.99$  $Date: 2011/11/15$

nobs = s.nobs;
nsys = s.nsys;
npar = s.npar;

if nargin<2
    i=nsys+npar;
end

if ischar(i)
    if ~isempty(s.names)
        [~,~,i] = intersect(i,s.names((nobs+1):end));
        if isempty(i)
            fprintf('Wrong variable name, choose among:\n');
            for k=nobs+(1:(nsys+npar-1))
                fprintf('%s, ',s.names{k});
            end
            fprintf('%s.\n',s.names{end});
            error('dfxiplot wrong variable name');
        end
    else
        error('names is not defined in CATS structure... can''t find appropriate variable.');
    end
elseif isfloat(i)
    if(i>(nsys+npar))||(i<1)
        error(['i should be between 1 and' num2str(nsys+napr)]);
    end
else
    error('i has to be a char string or an integer to reference an existing variable in CATS structure.');
end
clf;
barh(s.dfxs(i,1:nobs));
colormap(lines(nsys));
cmap = [.5*ones(1,3); circshift(lines(nsys),-1)];
hold on;
for j=1:nsys
    barh(nobs+j,s.dfxs(i,nobs+j),'FaceColor',cmap(j,:));
end
hold off;
box on;
axis([min(0,min(s.dfxs(i,:))*1.2) max(0,max(s.dfxs(i,:)))*1.2 0 length(s.dfxs(i,:))+1 ]);
grid off;
% set(gca,'Xtick',round(linspace(1,length(s.leverage),min(10,nobs+nsys))));
if isempty(s.names)
    set(gca,'YDir','Reverse','Ytick',1:(nobs+nsys),'YtickLabel',1:(nobs+nsys));
    grid on;
    set(gca,'YGrid','off');
else
    set(gca,'YDir','Reverse','Ytick',1:(nobs+nsys),'YtickLabel',s.names(1:(nobs+nsys)));
    PrettyFigureFormat;
    grid on;
    set(gca,'YGrid','off');
    %     rotateticklabel(gca,90);
    %     set(gca,'Visible','off');
end
if ~isempty(s.names)
    xlabel(sprintf('Standardized change in fitted value of %s',s.names{nobs+i}));
end

end