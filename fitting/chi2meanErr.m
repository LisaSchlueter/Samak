% chi2 function to estimate:
% - parMean (and it's uncertainty). This is e.g. the mean value of the endpoint, fitted individually by n runs
% input: 
% -parRun: values for parameter from the individuals fits (e.g. endpoint values)
% -CM: Covariance Matrix, which contains: different uncertainties on parameter and their correlations
function chi2 = chi2meanErr(parMean,parRunCM)
parRun = parRunCM{1};
CM = parRunCM{2};
chi2 = ((parRun-parMean)')*(CM\(parRun-parMean));
end
