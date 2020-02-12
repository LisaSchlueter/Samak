m     = 0.360 * 2 * 3600; % mean count rate per qU
qUmE0 = 1:.1:100;
runs  = 1;

% Poisson Simulation
p         = m + sqrt(m).*randn(runs,numel(qUmE0)); % raw=runs col=qU bins;
pstack    = sum(p,1);
pstackFit = fitdist(pstack','Normal');
fprintf(2,'Poissonian: mean=%g sigma=%g sigma/mean=%g percent (sqrt(m/n)=%g)\n',pstackFit.mu,pstackFit.sigma,pstackFit.sigma./pstackFit.mu*100,sqrt(mean(pstack,2)/runs));

x=2:0.5:10;
ratio=zeros(1,numel(x));
counter=0;
for i = 3
    counter=counter+1;
% Non-Poisson Simulation
alpha = i/100;
np         = m + alpha.*m.*randn(runs,numel(qUmE0)); % raw=runs col=qU bins;
npstack    = sum(np,1);
npstackFit = fitdist(npstack','Normal');
fprintf(2,'Non-Poissonian: mean=%g sigma=%g sigma/mean=%g percent\n',npstackFit.mu,npstackFit.sigma,npstackFit.sigma./npstackFit.mu*100);
% Ratio
ratio(counter)=npstackFit.sigma./pstackFit.sigma;
fprintf('Non-Poissonian (%g %%) / Poissonian Sigma = %g\n',i,ratio(counter));
end

%%
figure(1)
hns=histogram(npstack./(m*runs));
hold on
hs=histogram(pstack./(m*runs));
hs.BinWidth=hns.BinWidth;
hold off
PrettyFigureFormat
xlabel('Counts/qU');

%%
figure(2)
plot(x,ratio,'LineWidth',5,'Color',rgb('Amethyst'));
xlabel('Non-Poissonian Background in Percent (2h run)');
ylabel('\sigma_{np}/\sigma_{p}');
title('Non-Poissonian Background Correction Factor - KNM1 - 2h runs - 360 runs')
PrettyFigureFormat
