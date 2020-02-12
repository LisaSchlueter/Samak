function chi2 = Chi2E0_first(p,XA)
% Chi2 Function For Fit Tritium Beta Spectrum
% With Nuisance Parameters
% 
% Th. Lasserre, 2012-13
%
% Gaussian Xi2
% X = [ "column vector of bin centers"  | "column vector of data" | "column vector of uncertainties" | "... extra info" ]
% p = [ p1; p2; ... ] column vector of parameter to fit

X = XA{1};
A = XA{2};

A.ComputeTBDDS(...
    'mSq_bias',p(1),...
    'E0_bias',p(2),...
    'B_bias',p(3),...
    'N_bias',p(4)); 
A.ComputeTBDIS();

% Poisson deviance chi square.
chi2 = sum(((X(:,2) - A.TBDIS)./X(:,3)).^2);
            % counts            %uncertainties
end