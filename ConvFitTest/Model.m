function f = Model(p,x)
% Model function, p are model parameters, x are model inputs

% G. Mention - CEA Saclay, December 2012

Poisson = p(1).^x.*exp(1).^(-p(1))./gamma(x+1);
Gauss   = exp(-0.5.*((x-p(2))./p(3)).^2)./(p(3).*sqrt(2.*pi));

%f=p(1)*x+p(2);
%f = p(1)*(x+p(2)).^2;
%f = p(1) + p(2).*x + p(3).*x.^2 + p(4).*x.^3 ;
f = p(4).*conv(Poisson,Gauss,'same');

end