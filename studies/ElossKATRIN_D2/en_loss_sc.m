function f = en_loss_sc(x, amp1, pos1, sig1, amp2, pos2, sig2, amp3, pos3, sig3)
    % parametrized energy loss function, V. Hannen, March 2019, KATRIN CM36
    % Global fit of parametrized energy loss function to measured integral and t.o.f. data
    % Based on D2 scattering data
    
    T  = 18575;   % kinetic energy of incident electron [eV]
    Ei = 15.467;  % ionization energy for D2 [eV]
    
    f = (x <= 0).* 0  + ...
        (x >  0 & x <  Ei)  .*  (amp1 .* gaussian(x, pos1, sig1)   + amp2 .* gaussian(x, pos2, sig2)  + amp3 .* gaussian(x, pos3, sig3)) + ...
        (x >= Ei)  .* ((amp1 .* gaussian(Ei, pos1, sig1) + amp2 .* gaussian(Ei, pos2, sig2) + amp3 .* gaussian(Ei, pos3, sig3)) .* Dsigma(x, T, Ei) ./ Dsigma(Ei, T, Ei));
    f(isnan(f)) = 0;
end
