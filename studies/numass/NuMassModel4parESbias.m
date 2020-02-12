%
% Model function for spectrum minimization with minuit, 
% p are model fit parameters, com_opt are model input 
% values in pairs of "option" & "value".


%function f = Model(p,com_opt,DetModelType)
function f = NuMassModel4parESbias(p)

global A;

% Define Set of Bias to be applied
% A.mnuSq_Bias   = p(1);
% A.Q_Bias       = p(2);
% A.Bkg_Bias     = p(3);
% % norm
% A.mnu4Sq_Bias  = p(5);
% A.sin2T4_Bias  = p(6);

A.ComputeTBDDS(...
    'm_bias',p(1),...
    'Q_bias',p(2),...
    'B_bias',p(3),... 
    'TTES_bias',p(7)); 
A.ComputeTBDIS();
f = ((1+p(4)).*A.TBDIS);
end