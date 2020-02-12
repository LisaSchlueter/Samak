function chi2 = chi2lin(par,data)
model = par(1).*data(:,1)+par(2);
chi2 = sum((model-data(:,2)).^2./data(:,3).^2);
end