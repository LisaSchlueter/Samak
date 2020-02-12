function chi2 = KrLineChi2allPixels(p,X)
%
% Chi2 Function For Fit Krypton Spectrum
% All Pixels Defined in a Column Vector
% 
% Th. Lasserre, 2017
%
% Gaussian Xi2
% X = [ "column vector of bin centers"  | "column vector of data" | "column vector of uncertainties" | "... extra info" ]
% p = [ p1; p2; ... ] column vector of parameter to fit

nn = X(:,2)~=0; % indexes of non-null bins
x  = X(nn,1);    % bin centers of non-null bins
y  = X(nn,2);    % bin content of non-null bins
z  = X(nn,3);    % uncertainties bin content of non-null bins
m  =  KrLineModelallPixels(p);

% Gaussian
chi2 = sum(( y(1:numel(x)) - m(1:numel(x))).^2 ./ z(1:numel(x)).^2);
end