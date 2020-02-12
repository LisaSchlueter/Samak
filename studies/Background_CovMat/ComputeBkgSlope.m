function Bkg = ComputeBkgSlope(par,qU)
BkgFlat = par(1);
Slope   = par(2);
Bkg = arrayfun(@(x) Slope*x+BkgFlat,qU-18573.8);
end