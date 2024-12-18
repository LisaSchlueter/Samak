function deleteparplot(s)
%% Delete-1 parameters

XD=1:size(s.x_i,2);
YD=1:size(s.x_i,1);
diffxi = s.x_i-repmat(s.xhat,[1 size(s.x_i,2)]);
imagesc(XD,YD,diffxi./repmat(max(abs(diffxi),[],2),[1 size(s.x_i,2)]));
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
set(gca,'YTickLabel',{'E_0' 'B' 'N'})
title('Delete-1 parameters (centered and rescaled)');

end