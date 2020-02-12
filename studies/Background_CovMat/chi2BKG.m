function chi2 = chi2BKG(par,data)
qU = data(:,1);
BKG = ComputeBkgSlope(par,qU); % Computes array of background with slope
chi2 = sum((BKG-data(:,2)).^2./data(:,3).^2);
end