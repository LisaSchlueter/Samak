function chi2 = Chi2GaussMultiRing(p,X)
% Gaussian Xi2
% X = [ "column vector of bin centers"  | "column vector of data" | "column vector of uncertainties" | "... extra info" ]
% p = [ p1; p2; ... ] column vector of parameter to fit

% Th. Lasserre - CEA Saclay - December 2012

x1 = X(1:26,1);    % bin centers of non-null bins
y1 = X(1:26,2);    % bin content of non-null bins
z1 = X(1:26,3);    % uncertainties bin content of non-null bins

x2 = X(27:52,1);    % bin centers of non-null bins
y2 = X(27:52,2);    % bin content of non-null bins
z2 = X(27:52,3);    % uncertainties bin content of non-null bins

x3 = X(53:78,1);    % bin centers of non-null bins
y3 = X(53:78,2);    % bin content of non-null bins
z3 = X(53:78,3);    % uncertainties bin content of non-null bins

x4 = X(79:104,1);    % bin centers of non-null bins
y4 = X(79:104,2);    % bin content of non-null bins
z4 = X(79:104,3);    % uncertainties bin content of non-null bins

m1 = Model(p(1:2),x1);
m2 = Model(p(3:4),x2);
m3 = Model(p(5:6),x3);
m4 = Model(p(7:8),x4);

slopeError = 1e4;

chi2 = sum(( y1 - m1 ).^2 ./ z1.^2) + ...
       sum(( y2 - m2 ).^2 ./ z2.^2) + ...
       sum(( y3 - m3 ).^2 ./ z3.^2) + ...
       sum(( y4 - m4 ).^2 ./ z4.^2) + ...
       ((p(4)/p(3)-p(2)/p(1))/slopeError).^2 + ... 
       ((p(6)/p(5)-p(4)/p(3))/slopeError).^2 + ...
       ((p(8)/p(7)-p(6)/p(5))/slopeError).^2 + ...
       ((p(8)/p(7)-p(2)/p(1))/slopeError).^2;
end
