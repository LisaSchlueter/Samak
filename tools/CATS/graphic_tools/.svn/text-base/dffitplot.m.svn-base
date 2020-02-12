function dffitsplot(s)
%% Scaled change in fitted values. DFFITS (with studentized residuals)

bar(s.dffit(1:s.nobs));
colormap(lines(s.nsys+s.npar));
cmap = [.5*ones(1,3); circshift(lines(s.nsys+s.npar),-1)];
hold on;
for i=1:s.nsys
    bar(s.nobs+i,s.dffit(s.nobs+i),'FaceColor',cmap(i,:));
end
hold off;
box on;
axis([0 length(s.dffit)+1 1.2*max(abs(s.dffit))*[-1 1]]);
grid off;
set(gca,'Xtick',round(linspace(1,length(s.dffit),min(10,s.nobs+s.nsys))));
xlabel('Bin number');
ylabel('dffits');
title('Scaled change in fitted values. DFFITS');

end