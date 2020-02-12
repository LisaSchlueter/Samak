function chi2 = Chi2Gauss(p,X)
% Gaussian Xi2
% X = [ "column vector of bin centers"  | "column vector of data" | "column vector of uncertainties" | "... extra info" ]
% p = [ p1; p2; ... ] column vector of parameter to fit

% Th. Lasserre - CEA Saclay - December 2012

nn = X(:,2)~=0; % indexes of non-null bins
x = X(nn,1);    % bin centers of non-null bins
y = X(nn,2);    % bin content of non-null bins
z = X(nn,3);    % uncertainties bin content of non-null bins
% m = p(1)+p(2)*(x-p(3)).^2; % model
m = Model(p,x);

chi2 = sum(( y - m ).^2 ./ z.^2) ;

end
