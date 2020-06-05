%
% Relic Neutrino Spectrum Per Specie
% 
% Th. Lasserre 15/06/2014
%

% Call Class
clc ; clf ; clear all;  A = TritiumRelicNu('RNS_nTnu',1000);

% Relic Neutrino Spectrum
h0 = semilogx(A.RNS_Enu,A.RNS_Snu,'LineWidth',2,'Color','Black');
grid on;
xlabel('Energy (keV)','FontSize',14);
str = sprintf('dN/dE (per %.1s keV)',A.TnuStep);
ylabel(str,'FontSize',14);
title('Relic Neutrino Spectrum Per Specie','FontSize',14);
str = sprintf('Fermi-Dirac - 1.9 K - %.0f /cm^3 per specie',sum(A.RNS_Snu));
h = legend([h0],str);
%legend(h,'boxoff');
set(h,'FontSize',12);
PrettyFigureFormat;
publish_figure(1,'relicnuspectrum.eps');