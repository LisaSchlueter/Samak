function DeltaChi2 = GetDeltaChi2(cl,npar)
% Calculate delta chi2 for a : confiedence level (cl), number of parameters (npar)
% DeltaChi2 = GetDeltaChi2(cl,npar)
% Lisa, April 2020
if cl>1
    cl = cl/100;
end
DeltaChi2 = chi2inv(cl,npar); % inverse cumulative chi2 distribution
end