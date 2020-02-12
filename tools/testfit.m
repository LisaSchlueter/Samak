% Data
xdata = (-5:.1:5)';
S=1; A = 100; B = 500;
ydata_i  = A + B * 1/(2*pi*S^2).*exp(-xdata.^2/S^2);
ydatae = sqrt(ydata_i);

% Fit Function
myf = @(a,b,s,x) a + b .* (1/(2*pi*s^2).*exp(-x.^2/s^2));

% Store Fit Parameters
nfit = 100;
chi2v = zeros(1,nfit);
out = 'ON';

% Fit Options
%    'Robust','BiSquare',...
f = fittype(@(a,b,s,x) myf(a,b,s,x),...
    'independent', 'x');
% Fit Algo
for (k=1:1:nfit)
    ydata  = ydata_i + ydatae .* randn(numel(xdata),1);
    opts = fitoptions('StartPoint',[A,B,S],...
        'Method', 'NonlinearLeastSquares',...
        'Weights',1./ydata_i,...
        'Display','off');
    [fit1,gof,fitinfo] = fit(xdata,ydata,f,opts);
    chi2v(k)= gof.sse;
    switch out
        case 'ON'
            ci = confint(fit1,0.68);
            fprintf('A = %3g \t fit = %3g \t (%3g - 68%%) \n',A,fit1.a,ci(2,1)-ci(1,1));
            fprintf('B = %3g \t fit = %3g \t (%3g - 68%%) \n',B,fit1.b,ci(2,2)-ci(1,2));
            fprintf('S = %3g \t fit = %3g \t (%3g - 68%%) \n',S,fit1.s,ci(2,3)-ci(1,3));
            fprintf('Chi2 = %g \n',chi2v(k));
            fprintf(2,'---\n');
    end
end

% Plot
figure(1)
hdata = errorbar(xdata,ydata,ydatae,...
    'ks','MarkerSize',5,'MarkerFaceColor',.8*[1 1 1],'LineWidth',1);
errorbar_tick(hdata,200);
hold on;
plot(fit1,'r-',xdata,ydata,'k.')
hold off;
grid on
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
set(gca,'FontSize',12);
set(gca,'yscale','lin');
PrettyFigureFormat

figure(2)
nhist(chi2v);
PrettyFigureFormat