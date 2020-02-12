function chi2 = krl332minuit_chi2cov(p,allX)
%
% Chi2 Function For Fit Krypton Spectrum
% With Nuisance Parameters
% 
% Th. Lasserre, 2017
%
% Gaussian Xi2
% allX = {Data (Matrix), Object A}  -> Cell strucuture
% X = [ "column vector of bin centers"  | "column vector of data" | "column vector of uncertainties" | "... extra info" ]
% p = [ p1; p2; ... ] column vector of parameter to fit

%global FitCovErrMatrix;
X  = allX{1,1};  %Data Vector
nn = X(:,2)~=0;  % indexes of non-null bins
x  = X(nn,1);    % bin centers of non-null bins
y  = X(nn,2);    % bin content of non-null bins
z  = X(nn,3);    % uncertainties bin content of non-null bins
A  = allX{1,2};    %Krypton Object
FitCovErrMatrix = allX{1,3};
m  =  krl332minuit_modelint(p,A);

% Poisson deviance chi square.
% chi2 = -2*sum( y(1:numel(x)) - m(1:numel(x)) + y(1:numel(x)).*log(m(1:numel(x))./y(1:numel(x))) );                      

% Gaussian
chi2 = (y(1:numel(x)) - m(1:numel(x)))' * ( FitCovErrMatrix \ ( y(1:numel(x)) - m(1:numel(x))));

end