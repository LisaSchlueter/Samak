function chi2 = Chi2E0mpixCM(p,XACM)
% Chi2 Function For Fit Tritium Beta Spectrum
% With Nuisance Parameters
% 
% Th. Lasserre, 2012-13
%
% Gaussian Xi2
% X = [ "column vector of bin centers"  | "column vector of data" | "column vector of uncertainties" | "... extra info" ]
% p = [ p1; p2; ... ] column vector of parameter to fit

X = XACM{1};
Afit = XACM{2};
FCEM = XACM{3}; %FitCovErrMatrix

y = X(:,2);    % bin content of non-null bins

Afit.ComputeTBDDS(...
    'mSq_bias',p(1),...
    'E0_bias',p(2),...
    'N_bias',p(Afit.nPixels+3:Afit.nPixels+3+Afit.nPixels-1),...
    'B_bias',p(3:Afit.nPixels+2));
% p(3:A.nPixels+2)p(A.nPixels+3:A.nPixels+3+A.nPixels-1)
Afit.ComputeTBDIS();

m = Afit.TBDIS(:);

% Poisson deviance chi square.
% With pull term for the neutrino mass for first tritium studies
chi2 = (y - m)' * ( FCEM \ (y - m)) + (p(1)-0)^2/4;

end