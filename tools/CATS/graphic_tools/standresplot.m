function standresplot(s)
% Standardize residual plot of the fit performed by CATS
%
% resplot(s) plots the standardized residuals of the fit performed by CATS
% where s is the output result of CATS ( s = cats(...) ).
%
% Copyright 2005-2010 Guillaume MENTION, CEA Saclay
% $Revision: 0.97$  $Date: 2010/09/14$

nobs = length(s.yhat);
nsys = 0;
npar = length(s.xhat);

clf;
bar(s.standres,'facecolor',rgb('CadetBlue'));
title('Residuals')
ylabel('R');


end