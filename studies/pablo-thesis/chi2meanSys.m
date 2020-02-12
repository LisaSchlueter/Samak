function chi2 = chi2meanSys(parMean,parRunCM)
parRun = parRunCM{1};
CM = parRunCM{2};
chi2 = ((parRun-parMean)')*(CM\(parRun-parMean));
a=1;
end