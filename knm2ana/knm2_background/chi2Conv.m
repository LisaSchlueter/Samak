function chi2 = chi2Conv(par,data)

% par 1 == sigma poisson
% par 2 == sigma gauss
% par 3 == mu gauss
% data(:,3) == counts occurence

xConv       = linspace(-max(data(:,1))*1.2,max(data(:,1)*1.2,1e3));
GaussDist   = gaussian(xConv,par(3),par(2))./sqrt(2*pi*par(2)^2);
PoissDist   = gaussian(xConv,par(1),par(1))./sqrt(2*pi*par(1)^2);
ConvDist    = conv(GaussDist,PoissDist); 

chi2 = sum((model-data(:,2)).^2./data(:,3).^2);
end