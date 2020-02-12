% Model function for spectrum minimization with minuit, 
% p are model fit parameters, com_opt are model input 
% values in pairs of "option" & "value".

function f = model_tbdallpixelssim(p,A)

A.ComputeTBDDSallPixels(...
    'mSq_bias',p(1),...
    'E0_bias',p(2),...
    'Norm_BiasPixels',p(3:A.nPixels+2));
A.ComputeTBDISallPixels();
f = (reshape(A.TBDISallPixels(:,1:A.nPixels),A.nPixels*A.nqU,1));
end