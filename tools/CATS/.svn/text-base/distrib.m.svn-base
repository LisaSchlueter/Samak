function d = distrib(n,nobs)

x=linspace(10/(2*nobs),10-10/(2*nobs),nobs)';
sigma_b2b = 0.01;
t_acc = .10;
s_acc = .10;
t_li9 = .03;
s_li9 = .5;
t_prc = .10;
s_prc = .3;

norm = 15000;

nu = (x>1).*(x-1).^2.*exp(-(x-1));
ampb = (x>1).*(sin(1.27*2.5e-3*1.05e3./(x+0.782)).^2);
nu = nu/sum(nu);
acc = (x>.5).*exp(-x);
acc = acc/sum(acc);
li9 = (x>1).*(x<10).*(x-1).*(10-x);
li9 = li9/sum(li9);
prc = (1+.2*randn(size(x))).*(x>.5).*ones(size(x));
prc = prc/sum(prc);
s2t13 = .1;

d = zeros(n,1);

for i=1:n
y=norm*(nu.*(1-s2t13*ampb).*(1+1./(sqrt(norm*nu.*(1-s2t13*ampb))+eps).*randn(size(x))+sigma_b2b*randn(size(x)))) + ...
    norm*(t_acc*acc.*(1+s_acc*randn(size(x))) + ...
    t_li9*li9.*(1+s_li9*randn(size(x))) + ...
    t_prc*prc.*(1+s_prc*randn(size(x))));

A=-norm*ampb.*nu;
S=norm*[nu acc li9 prc];
Winv = diag(1./[.03 s_acc s_li9 s_prc].^2);
x0=[1 t_acc t_li9 t_prc]';
%List = [1 2]';
Vinv=diag(1./(abs(y)+(norm*nu.*(1-s2t13*ampb)*sigma_b2b).^2+eps));

nsys = size(S,2);
npar = size(A,2);

s=cats(y,A,'Vinv',Vinv,'S',S,'Winv',Winv,'x0',x0);

d(i) = s.xhat(5);

end
