% Model function for spectrum minimization with minuit, 
% p are model fit parameters, com_opt are model input 
% values in pairs of "option" & "value".

function f = model_hvdata_pixeldiff(p,A)

%global A;
A.ComputeKrDS(...
    'E_bias',p(1),...
    'W_bias',p(2),...
    'Phi0_bias',p(3)); 
f = A.KrDS;
end