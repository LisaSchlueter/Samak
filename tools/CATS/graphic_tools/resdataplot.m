function resdataplot(s)
hold on;
if(RainbowFlag)
    cmap = jet(nobs);
    resid = s.r(1:nobs);
    for i=1:nobs,
        plot(y(i),resid(i),'ko','LineWidth',1,'MarkerFacecolor',cmap(i,:),...
        'MarkerSize',8);
    end
else
    plot(y,s.r(1:nobs),'bo');
end
hold off;
box on;
% Think more about this plot...
%hold on; plot(S*x0,s.r(1:nobs),'rx'); hold off;
grid on;
xlabel('N_{evt}');
ylabel('R(y)');
title('Residuals vs Data');
end
