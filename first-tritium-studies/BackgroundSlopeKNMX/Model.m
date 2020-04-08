function f = Model(p,x)
% Model function, p are model parameters, x are model inputs
f = p(1) + (x-18.5737) .* p(2);
end