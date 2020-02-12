function chi2 = NuMassChi2E0_Cov(p,XA)
% Chi2 Function For Fit Tritium Beta Spectrum
% With Nuisance Parameters
% 
% Th. Lasserre, 2012-13
%  L.Schlueter, 05/2018
% Gaussian Xi2
% X = [ "column vector of bin centers"  | "column vector of data" | "column vector of uncertainties" | "... extra info" ]
% p = [ p1; p2; ... ] column vector of parameter to fit

% CovMat

X = XA{1}; % Data
A = XA{2}; % Modell Object
FitCovMatFrac = XA{3}; % Fractional Covariance Matrix

nn = X(:,2)~=0; % indexes of non-null bins
x = X(nn,1);    % bin centers of non-null bins
y = X(nn,2);    % bin content of non-null bins
z = X(nn,3);    % uncertainties bin content of non-null bins

A.ComputeTBDDS(...
    'mSq_bias',p(1),...
    'E0_bias',p(2),...
    'B_bias',p(3),...
    'N_bias',p(4)); 
A.ComputeTBDIS('IStype','SIMP');
m = A.TBDIS;
FitCovMat = FitCovMatFrac.*m.*m';

%chi square
chi2 = (y-m)' * (FitCovMat \ (y-m)) + (p(1)-0)^2/4;

end