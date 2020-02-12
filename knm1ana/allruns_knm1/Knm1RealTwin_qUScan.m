function [RT_parqU, RT_errqU, RT_chi2qU, RT_dofqU] = Knm1RealTwin_qUScan(Real,Twin)
% Stacked Fit - Neutrino Mass Fixed

Real.fixPar = '1 5 6 7 8 9 10 11';
[RT_parqU.R, RT_errqU.R, RT_chi2qU.R, RT_dofqU.R] = Real.qUScan('saveplot','ON');

Twin.fixPar = '1 5 6 7 8 9 10 11';
[RT_parqU.T, RT_errqU.T, RT_chi2qU.T, RT_dofqU.T] = Twin.qUScan('saveplot','ON');

end

