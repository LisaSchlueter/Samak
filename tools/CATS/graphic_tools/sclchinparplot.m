function sclchinparplot(s)
%% Scale change in parameters

XD=1:size(s.dfxs,2);
YD=1:size(s.dfxs,1);
imagesc(XD,YD,s.dfxs);
set(gca,'Xtick',round(linspace(1,length(XD),min(8,s.nobs+s.nsys))));
set(gca,'YTick',round(linspace(1,length(YD),min(8,s.nsys+s.npar))));
caxis([-1 1]);
cm=0:1/9:1;
my_cm_blue=[1-cm;1-cm;ones(size(cm))]';
my_cm_red=[ones(size(cm));1-cm;1-cm]';
colormap([flipud(my_cm_blue);[1 1 1];my_cm_red]);
colorbar;
set(gcf,'Color',[1 1 1]);
colorbar;
set(gca,'XAxisLocation','top');
grid on;
%set(gca,'YTickLabel',{'nu' 'acc' 'li9' 'prc' 'dm2' 's2t13'})
title('Scale change in parameters');

end