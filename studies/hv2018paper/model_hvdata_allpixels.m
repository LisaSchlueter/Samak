% Model function for spectrum minimization with minuit, 
% p are model fit parameters, com_opt are model input 
% values in pairs of "option" & "value".

function f = model_hvdata_allpixels(p,A)

%global A;

A.ComputeKrDSallPixels(...
    'E_Bias',p(1),...
    'W_Bias',p(2),...
    'Phi0_Bias', p(3:A.nPixels+2),...
    'B_Bias', p(A.nPixels+3:2*A.nPixels+2));

A.ComputeKrISallPixels();
f = reshape(A.KrISallPixels,A.nPixels*A.nqU,1);
end