function dsdW = Dsigma(E, T, B)
% Unnormalized single differential ionization cross section.
% E = electron energy loss [eV]
% Ionization energy B for H2, D2 and T2 = 15.426 eV, 15.467 eV, 15.487 eV
% Minimal energy loss: Elossmin = B
% Maximal energy loss: Elossmax = (T - B)/2

% Binary encounter dipole model, eq. 44 in Kim, Rudd, Phys. Rev. A 50, 3954 (1994).
W  = E - B;   % secondary electron energy [eV]
w  = W./B;
t  = T./B;
y  = 1./(w+1); % = B/E for expansion of oscillator strength
N  = 2.;      % number of electrons in H2 molecule
Ni = 1.173;   % eq. 39 and table I.

dsdW = (Ni./N - 2)./(t+1) .* ( y + 1./(t-w) ) + (2 - Ni./N) .* (y.^2 + 1./(t-w).^2) + log(t).*y./N .* dfdwH2(y);

end
