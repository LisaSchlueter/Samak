%        par[0] amplitude of gaussian excitation peak
par0   =    +1.84234012E+03;  % error 8.25624701E+00
%        par[1] position of gaussian excitation peak
par1   =    +1.27186562E+01;  % error +2.54038843E-03
%        par[2] sigma of gaussian excitation peak
par2   =    +5.61673258E-01;  % error +1.87857676E-03
%        par[3] transition point between gaussian peak and line
par3   =    +1.36732536E+01;  % error +6.64567708E-03
%        par[4] transition point between line and tail
par4   =    +1.48591947E+01;  % error +4.86608385E-02
%        par[5] amplitude of tail 1/(DeltaE)**2
par5   =    -3.14460383E+03;  % error +1.21982053E+03
%        par[6] amplitude of tail 1/(DeltaE)**3
par6   =    +1.69971770E+06;  % error +2.30037595E+04
%        par[7] amplitude of twofold scattering
par7   = 0;

par = [par0 par1 par2 par3 par4 par5 par6 par7]; l1=par(4); l2=par(5);
gauss  = @(x,pos,sigma) exp(-(x-pos).*(x-pos)./(2.*sigma.*sigma));
tail   = @(x,square,cube) 0 + (x>0).*(square./(x+1e-99).^2+cube./(x+1e-99).^3);
y1     = par(1).*gauss(l1, par(2), par(3));
y2     = tail(l2,par(6),par(7));

f1scat = @(x,par) ...
    (x<l1).*par(1).*gauss(x,par(2),par(3)) + ...
    (x>=l1).*( (x<=l2).*(y1.*(l2-x) + y2.*(x-l1)./(l2-l1)) + (x>l2).*tail(x,par(6),par(7)));
f1scatn = @(e) f1scat(e,par)/simpsons(e,f1scat(e,par)); %normalized
x=-100:0.1:100;
p=Plot(x,f1scatn(x));
p.XLabel='Electron Kinetic Energy (eV)';
p.XLabel='pdf';
p.Title='KATRIN E-Loss Function - gaussian + single line with tail';
PrettyFigureFormat