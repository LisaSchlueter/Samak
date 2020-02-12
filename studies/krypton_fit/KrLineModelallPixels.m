% Model function for spectrum minimization with minuit, 
% p are model fit parameters, com_opt are model input 
% values in pairs of "option" & "value".

function f = KrLineModelallPixels(p)

global A;

A.ComputeKrDSallPixels(...
    'E_bias',p(1),...
    'W_bias',p(2));
% ...
%     'L3_32_Phi0allPixels',[100*ones(1,148)],...
%     'K32_Phi0allPixels',[100*ones(1,148)]);

A.ComputeKrISallPixels();
%A.ComputeKrISallPixels('BKG_RateSecallPixels',[1*ones(1,148)]);
f = (1+p(3))*cell2mat(reshape(A.KrISallPixels,A.nPixels*A.nqU,1));
end