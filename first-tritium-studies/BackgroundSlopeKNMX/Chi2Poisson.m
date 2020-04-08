
function chi2 = Chi2Poisson(p,X)
% Poisson deviance function: (asymptotically close to chi2 distribution)
% Chi2Poisson = -2*sum( y-m(p) + y*log(m(p)/y) )
% X = [ "column vector of bin centers"  | "column vector of data" | "column vector of uncertainties" | "... extra info" ]
% p = [ p1; p2; ... ] column vector of parameter to fit

% G. Mention - CEA Saclay - December 2012

nn = X(:,2)~=0; % indexes of non-null bins
x = X(nn,1);    % bin centers of non-null bins
y = X(nn,2);    % bin content of non-null bins
% m = p(1)+p(2)*(x-p(3)).^2; % model
m = Model(p,x);

chi2 = -2*sum( y - m + y.*log(m./y) ); % Poisson deviance chi square.

end
