function chi2 = Chi2GaussKNM12(p,X)
% Gaussian Xi2
% X = [ "column vector of bin centers"  | "column vector of data" | "column vector of uncertainties" | "... extra info" ]
% p = [ p1; p2; ... ] column vector of parameter to fit

% Th. Lasserre - CEA Saclay - December 2012

%knm1
x1 = X(1:7,1);    % bin centers of non-null bins
y1 = X(1:7,2);    % bin content of non-null bins
z1 = X(1:7,3);    % uncertainties bin content of non-null bins

%knm2
%knm2b = 8:numel(X(:,1));
knm2b = [8 14];
x2 = X(knm2b,1);    % bin centers of non-null bins
y2 = X(knm2b,2);    % bin content of non-null bins
z2 = X(knm2b,3);    % uncertainties bin content of non-null bins

m1 = Model(p(1:2),x1);
m2 = Model(p(3:4),x2);

slopeError = 1e-3;

chi2 = sum(( y1 - m1 ).^2 ./ z1.^2) + ...
       sum(( y2 - m2 ).^2 ./ z2.^2) + ...
       ((p(4)/p(3)-p(2)/p(1))/slopeError).^2;
end
