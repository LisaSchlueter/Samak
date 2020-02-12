% Model function for spectrum minimization with minuit, 
% p are model fit parameters, com_opt are model input 
% values in pairs of "option" & "value".

function f = tbd_modelint(p,A)
%A.qu = data 
%A.computephasespace

A.ComputeTBDDS(...
    'mSq_bias',p(1),...
    'E0_bias',p(2),...
    'B_bias',p(3),...
    'N_bias',p(4)); 
A.ComputeTBDIS('IStype','TRAPEZ');
f = A.TBDIS;
end
