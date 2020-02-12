function plotcorh(s)
% Correlation matrix plot of the fit performed by CATS
%
% plotcorh(s) plots the correlation matrix of the fit performed by CATS
% where s is the output result of CATS ( s = cats(...) ).
%
% Copyright 2005-2010 Guillaume MENTION, CEA Saclay
% $Revision: 0.97$  $Date: 2010/09/14$

clf;
nD=length(s.hatmat);
XD=1:nD;
YD=1:nD;
imagesc(XD,YD,sqrt(diag(1./diag(s.hatmat)))*s.hatmat*sqrt(diag(1./diag(s.hatmat))));
set(gca,'XTick',XD)
set(gca,'XTickLabel',{num2str(XD')})
set(gca,'YTick',YD)
set(gca,'YTickLabel',{num2str(YD')})
%set(gca,'YLabel',YD,YD-.5);
caxis([-1 1]);
cm=0:1/9:1;
my_cm_blue=[1-cm;1-cm;ones(size(cm))]';
my_cm_red=[ones(size(cm));1-cm;1-cm]';
colormap([flipud(my_cm_blue);[1 1 1];my_cm_red]);
colorbar;
set(gcf,'Color',[1 1 1]);
% set(gca,'GridLineStyle','-'); set(gca,'MinorGridLineStyle','-');
% set(gca,'XGrid','on'); set(gca,'XMinorGrid','on'); set(gca,'YGrid','on');
% set(gca,'YMinorGrid','on');
set(gca,'XAxisLocation','top');