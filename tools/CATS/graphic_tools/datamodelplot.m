function datamodelplot(s)
% to debug

grid on;
subplot(1,5,[1 4]);
dataplot = [s.S(1:s.nobs,5)*s.xhat(5) s.S(1:s.nobs,4)*s.xhat(4) s.S(1:s.nobs,3)*s.xhat(3)...
    s.S(1:s.nobs,2)*s.xhat(2) s.S(1:s.nobs,1)*s.xhat(1)+A(1:s.nobs,1)*s.xhat(s.nsys+1)];
bar(x,dataplot,'stacked');
colormap(flipud(lines(5)));
hold on;
stairs([x; x(end)+(x(end)-x(end-1))]-10/(2*s.nobs),...
    [s.yhat(1:s.nobs); s.yhat(s.nobs)],'Color',.5*ones(3,1),'LineWidth',2);
stairs([x; x(end)+(x(end)-x(end-1))]-10/(2*s.nobs),...
    [s.yhat(1:s.nobs)-A(1:s.nobs,1)*s.xhat(s.nsys+1);...
    s.yhat(s.nobs)-A(s.nobs,1)*s.xhat(s.nsys+1)],'Color','k','LineWidth',2);
if(RainbowFlag)
    cmap = jet(s.nobs);
    for i=1:s.nobs,
        errorbar(x(i),y(i),1/sqrt(Vinv(i,i)),'ks','LineWidth',1,'MarkerFacecolor',cmap(i,:),...
        'MarkerSize',4);
    end
else
    errorbar(x,y,1./sqrt(diag(Vinv)),'ks','LineWidth',1,'MarkerFacecolor','y',...
        'MarkerSize',4);
end
hold off;
axis([0 10 min(min(min(dataplot,[],1))*1.2,0) max(y*1.2)]);
set(gca,'Xtick',0:10);
xlabel('E_{vis} in MeV');
ylabel('N_{evt}');
title('Data and Model');
legend('Dm2 Systematic','Proton Recoils','Li9','Accidentals','\nu signal',...
    'Total Oscillation Model','No oscillation Model','Data','Location',...
    'EastOutside');

subplot(1,5,5);
cmap = [.5*ones(1,3); circshift(lines(5),-1)];
hold on;
ini_p = [x0;s2t13]; % to check behaviour only with s2t13. In principle,
% no prior on its value...
pulls = zeros(s.nsys+s.npar,1);
for i=1:length(s.xhat)
    pulls(i) = (s.xhat(i)-ini_p(i))./sqrt(s.covx(i,i));
    bar(i,pulls(i),'FaceColor',cmap(i,:));
end
%errorbar(1:length(x0),x0,1./sqrt(diag(Winv)),'ko','LineWidth',2,'MarkerFac
%ecolor','y','MarkerSize',10);
axis([0 length(x0)+2 -1.2*max(abs(pulls)) 1.2*max(abs(pulls))]);
box on;
hold off;
title('Reduced Pulls');


end