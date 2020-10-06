function chi2 = chi2linMean(par,data)
model = par(1).*data(:,1)+wmean(data(:,2),1./data(:,3).^2);
chi2 = sum((model-data(:,2)).^2./data(:,3).^2);
end