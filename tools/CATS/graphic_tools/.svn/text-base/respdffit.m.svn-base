function [h1 h2] = respdffit(s,alpha,bw)
% [h1 h2] = respdffit(s) plots the pdf of the standardized residuals of the
% fit performed by CATS and fit it with a normal distribution.
%  
%
if nargin<2
    alpha = erfc(1/sqrt(2));
end
if nargin<3
    bw=1;
end

clf;

[hh xx]=hist(s.standres,-5:bw:5);
hb = bar(xx,hh/sum(hh),'BarWidth',1,'FaceColor','b');

hold on;
x=linspace(-5,5,1000);
[muhat,sigmahat,muci,sigmaci] = normfit(s.standres,alpha);
hp = plot(x,normpdf(x,muhat,sigmahat)*bw,'r','LineWidth',2);
hold off;
PrettyFigureFormat;
set(gca,'Xtick',-5:1:5);
% grid on;
xlabel('Standardized residuals pdf');
text(-5,.4*bw,sprintf('\\mu = %1.2f ± %1.2f\n\\sigma = %1.2f ± %1.2f\n',...
    muhat,diff(muci)/2,sigmahat,diff(sigmaci)/2),...
    'BackgroundColor','none','FontSize',14,...
    'FontName','Helvetica','FontWeight','bold');
axis([-6 6 0 .6*bw]);
switch nargout
    case 1
        h1 = hb;
    case 2
        h1 = hb;
        h2 = hp;
end

end