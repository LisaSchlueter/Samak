function Knm1RealTwin_FitStack(Real,Twin)
% Stacked Fit - Neutrino Mass Fixed

Real.fixPar = '1 5 6 7 8 9 10 11';
Real.Fit; Real.PlotFit('Mode','Rate','saveplot','png');

Twin.fixPar = '1 5 6 7 8 9 10 11';
Twin.Fit; Twin.PlotFit('Mode','Rate','saveplot','png');

end
