function f = dfdwH2(y)
% Differential dipole oscillator strength for H2; table I in  Kim, Rudd, Phys. Rev. A 50, 3954 (1994).

c = 1.1262;
d = 6.3982;
e = -7.8055;
f = 2.1440;

f =  c .* y.^3 + d .* y.^4 + e .* y.^5 + f .* y.^6;

end