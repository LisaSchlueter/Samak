%
% Model of a Kr83m Single Line
% Model function for spectrum minimization with minuit 
% Thierry Lasserre, Last Updated: Jan. 2018
% 
% E:    Line Position  - Lorentian
% W:    Line Width     - Lorentian
% Phi0: Line Amplitude - cps
% B:    Offset         - cps (integral spectrum)
% 
function f = KrLineModel4par(p,A)
A.ComputeKrDS(...
    'E_bias',p(1),...
    'W_bias',p(2),...
    'Phi0_bias',p(3),...
    'B_bias',p(4));
A.ComputeKrIS();
f = ((1+p(5)).*A.KrIS);
end