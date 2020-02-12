function chi2 = NuMassChi2E0_CovB(p,XA)
% Chi2 Function For Fit Tritium Beta Spectrum
% With Nuisance Parameters
% 
% Th. Lasserre, 2012-13
%
% Gaussian Xi2
% X = [ "column vector of bin centers"  | "column vector of data" | "column vector of uncertainties" | "... extra info" ]
% p = [ p1; p2; ... ] column vector of parameter to fit

X = XA{1};         % Data
A = XA{2};         % Model (Object)
FitCovMat = XA{3}; %Covariance Matrix
FitLB = XA{4};     % Lower  Bin for Fit
FitHB = XA{5};     % Higher Bin for Fit
FitCovMatBfrac = XA{6}; %Covariance Matrix

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

% Background + Signal Covariance Matrix
ScaleB      = A.BKG_RateSec*A.qUfrac.*A.TimeSec;
FitCovMatTB = FitCovMat(FitLB:end,FitLB:end) + ( ScaleB(FitLB:end).*FitCovMatBfrac(FitLB:end,FitLB:end).*ScaleB(FitLB:end)');

% chi square.
chi2 = (y(FitLB:numel(x)) - m(FitLB:numel(x)))' * ( FitCovMatTB \ ( y(FitLB:numel(x)) - m(FitLB:numel(x)))) + (p(1)-0)^2/4;% + (p(3)-0)^2/0.005^2;


end