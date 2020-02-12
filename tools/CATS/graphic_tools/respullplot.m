function respullplot(s)
% to debug - sept 2010

colormap(lines);
set(gcf,'OuterPosition',[XY_scr large_figure_size]);
subplot(1,3,1);
bar(s.r(1:nobs));
axis([0 length(s.r(1:nobs))+1 min(s.r(1:nobs))*1.2 max(s.r(1:nobs))*1.2]);
grid off;
set(gca,'Xtick',round(linspace(1,length(s.r(1:nobs)),10)));
xlabel('Bin number');
ylabel('R');
title('Residuals');

subplot(1,3,2);
qqplot(s.r(1:nobs));
box on;
grid on;
title('Residual QQ-plot');

subplot(1,3,3);
colormap(lines(5));
cmap = [.5*ones(1,3); circshift(lines(5),-1)];
hold on;
for i=1:nsys
    bar(i,s.r(nobs+i),'FaceColor',cmap(i,:));
end
hold off;
box on;
axis([0 length(s.r(nobs+(1:nsys)))+1 ...
    min(s.r(nobs+(1:nsys)))*1.2 max(s.r(nobs+(1:nsys)))*1.2]);
grid off;
set(gca,'Xtick',round(linspace(1,length(s.r(nobs+(1:nsys))),min(8,nsys))));
xlabel('Bin number');
ylabel('R');
title('Pulls');

end