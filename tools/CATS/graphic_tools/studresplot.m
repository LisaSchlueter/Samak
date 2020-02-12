function studresplot(s)
% Studentized residual plot of the fit performed by CATS
%
% resplot(s) plots the studentized residuals of the fit performed by CATS
% where s is the output result of CATS ( s = cats(...) ).
%
% Copyright 2005-2010 Guillaume MENTION, CEA Saclay
% $Revision: 0.97$  $Date: 2010/09/14$

nobs = length(s.yhat);
nsys = 0;
npar = length(s.xhat);

clf;
bar(s.studres,'facecolor',rgb('CadetBlue'));
title('Residuals')
ylabel('R');


end