function chi2 = tbd_chi2Poisson(p,XA)
%
% Chi2 Function For Fit Krypton Spectrum
% With Nuisance Parameters
% 
% Th. Lasserre, 2017
%
% Gaussian Xi2
% X = [ "column vector of bin centers"  | "column vector of data" | "column vector of uncertainties" | "... extra info" ]
% p = [ p1; p2; ... ] column vector of parameter to fit

% New implementation without global variable
X = XA{1};
A = XA{2};

nn = X(:,2)>=0;  % indexes of non-null bins
x  = X(nn,1);    % bin centers of non-null bins
y  = X(nn,2);    % bin content of non-null bins
z  = X(nn,3);    % uncertainties bin content of non-null bins
m  =  tbd_modelint(p,A); m(isnan(m))=1e-9; y(y==0)=1e-9;
% Poisson deviance chi square.

 chi2 = -2*sum( y(1:numel(x)) - m(1:numel(x)) + y(1:numel(x)).*log(m(1:numel(x))./y(1:numel(x))) );     
 chi2 = chi2 + (p(1).^2)/5^2; % 5eV
 chi2 = chi2 + (p(2).^2)/2^2;  % 5eV
 chi2 = chi2 + (p(3).^2)/2^2; % 1cps
 chi2 = chi2 + (p(4).^2)/2^2; % 100%
end