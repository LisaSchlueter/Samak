
function f = KrLineModel7par_sat(p)

global A;

A.ComputeKrDS(...
    'E_Bias',p(1),...
    'W_Bias',p(2),...
    'Phi0_Bias',p(3),...
    'B_bias',p(4),...
    'E_Bias_s3', p(6),...
    'Phi0_s3_Bias', p(7));
A.ComputeKrIS();
f = ((1+p(5)).*A.KrIS);
end