% Model function for spectrum minimization with minuit, 
% p are model fit parameters, com_opt are model input 
% values in pairs of "option" & "value".

function f = KrLineModelDiff4par(p)

global A;

A.ComputeKrDS(...
    'E_bias',p(1),...
    'W_bias',p(2),...
    'B_bias',p(4)); 
f = ((1+p(3)).*A.KrDS);
end