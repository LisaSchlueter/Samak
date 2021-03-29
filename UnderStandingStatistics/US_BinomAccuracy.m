% understanding statistics better :)
mu    = 0;             
sigma = 1;
nSamples   = 500;   
nReSamples = 1e4;

x = mu+sigma.*randn(nSamples,nReSamples); % sampling: repeat experiment nReSamples times (1 experiment: nSamples observations of x) 

%% 1) sanity check: known case of sample mean and error of the mean
% because of central limit theorem sample mean is gaussian distributed
% ErrorOfMean = std(mean(x));
% fprintf('Expected error of the mean = %.4f \n',sigma/sqrt(nSamples));
% fprintf('MC error of the mean       = %.4f \n',ErrorOfMean);

%% 2) fraction that is smaller than x_bf
% number of x<x_bf follows not gaussian, but binomial dsitribution (discrete!)
% fraction of x<x_bf follows which distribution? Something like relative binomial...test this
% 
x_bf = -1; % threshold x
Frac_bf = sum(x<=x_bf)./nSamples; % fraction of x smaller than x_bf
Num_bf = sum(x<=x_bf);             % number of x smaller than x_bf
ErrorOfFrac = std(Frac_bf);        % sample error of fram

% estimated from binomial distrbution
p = normcdf(x_bf,mu,sigma); % known probability for x<=x_bf as input for binomial dist.
mu_bin     = nSamples*p;
var_bin    = nSamples*p.*(1-p);
relSigma_bin = sqrt(var_bin./nSamples^2);
fprintf('-----------------------------%.0f samples -------------------------------------------------\n',nSamples)
fprintf('Expected numbers of x<=x_bf from binomial  = %.1f +- %.1f (rel. err %.1f %%) \n',mu_bin,sqrt(var_bin),1e2*relSigma_bin);
fprintf('Sample MC numbers of x<=x_bf from binomial = %.1f +- %.1f (rel. err %.1f %%) \n',mean(Num_bf),std(Num_bf),1e2*std(Num_bf)/nSamples);
fprintf('----------------------------------------------------------------------\n')
fprintf('Expected fraction (binomial distribution) %.1f%% +- %.1f%%\n',100*mu_bin/nSamples,100*sqrt(var_bin)/nSamples);
fprintf('MC fraction:                              %.1f%% +- %.1f%%\n',100*mean(Frac_bf),100*ErrorOfFrac);
fprintf('--------------------------------------------------------------------------------------------\n')
