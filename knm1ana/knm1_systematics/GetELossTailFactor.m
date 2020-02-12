function out = GetELossTailFactor(obj,Emax)

% compute energy loss function: first scattering for full range
maxE_full = 9288;
minE_full=-maxE_full; NbinE_full = (maxE_full-minE_full)/0.1;
E_full = linspace(minE_full,maxE_full,NbinE_full);

[~, Eloss_full] = obj.ComputeELossFunction('E',E_full);

%integrate until Emax
[~, index_up]  = min(abs(E_full-Emax));
[~, index_low] = min(abs(E_full+Emax));

Eloss_full1scat = Eloss_full(1,:); %first scattering
out = simpsons(E_full(index_low:index_up),Eloss_full1scat(index_low:index_up));
end