% Convert Delta chi^2 from KSN-2 into confidence levels or sigmas


%DeltaChi2 = ConvertCLStd('Mode','CL2Sigma','nPar',2,'CL',0.95);
DeltaChi2 =0.74;

Sigma1 = @(DeltaChi2) sqrt(-2.*log(1-chi2cdf(DeltaChi2,2)));
Sigma2 = @(p) sqrt(-2.*log(1-p));

Sigma1(DeltaChi2)