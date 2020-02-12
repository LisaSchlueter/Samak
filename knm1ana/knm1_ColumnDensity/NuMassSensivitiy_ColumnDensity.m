A = MultiRunAnalysis('RunList','KNM1','chi2','chi2CM','exclDataStart',17); %take 17 for covariance matrix normalization

A.ModelObj.ComputeTBDDS; A.ModelObj.ComputeTBDIS;
A.RunData.TBDIS = A.ModelObj.TBDIS;
A.RunData.TBDISE = sqrt(A.ModelObj.TBDIS);

%%
A.fixPar = ' 5 6 7 8 9 10 11';
A.exclDataStart = 17; %problem: for lower ranges uncertainty estimation seems to fail, due to negativ nu-mass??

A.chi2 = 'chi2CM';
A.ComputeCM;
A.Fit;
mNuSqCM = A.FitResult.par(1);
mNuSqErrCM = A.FitResult.err(1);


A.chi2 = 'chi2Stat';
A.Fit;
mNuSqStat = A.FitResult.par(1);
mNuSqErrStat = A.FitResult.err(1);
%%
fprintf('%.0feV range \n',round(A.ModelObj.qU(A.exclDataStart)-18575));
fprintf('sensitivity mNu 90%%C.L. = %.3f (stat) \n',sqrt(mNuSqErrStat*1.64));
fprintf('sensitivity mNu 90%%C.L. = %.3f (stat+sys) \n',sqrt(mNuSqErrCM*1.64));


