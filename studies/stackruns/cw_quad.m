x = 20:10:200;
f = @(qu) 3000./(.1*(qu+10))+10;
y = f(x) + randn(1,numel(x)).*sqrt(f(x));
xmean = mean(x);
ymean = mean(y);
fmean  = f(mean(x));

% Solve
g = @(qu) f(qu) - ymean;
sol = fsolve(g,xmean);
disp(sol)

figure(1)
hf = plot(x,f(x));
hold on
hy=plot(x,y,'LineStyle','--','Marker','s','Color','Green');
hfmean = plot(xmean,fmean,'d','Color','red','MarkerSize',10);
hymean = plot(xmean,ymean,'d','Color','blue','MarkerSize',10);
hsol = plot(sol,ymean,'d','Color','black','MarkerSize',10);
hold off
xlabel('qu');
ylabel('N');
legend([hf hy hfmean hymean hsol],'f','data (qu,N)','f(mean data qu)','mean data N','solve N(qu stack)=mean data N');
PrettyFigureFormat;



