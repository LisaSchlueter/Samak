function [par, err, chi2min,dof,errmat] = linFit(x,y,yErr)
% input: x,y,sigma(y)
% output: 1. coefficients and 2. uncertainties on coefficients of a linear fit,
%         3. goodness of fit (chi2), 4. degrees of freedom
tmparg = sprintf(['set pri 1; min ; minos '],'');
Data = [x,y,yErr];
parInit = [0 0];
Args   = {parInit, Data, '-c', tmparg};
[par, err, chi2min, errmat] = fminuit('chi2lin',Args{:});
dof = numel(y)-numel(parInit);
end

