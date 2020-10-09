%
% Relic Neutrino Capture on Tritium
% Th. Lasserre - June 2014
%
addpath('../nukev')

% Call Class
clc ; clf ; clear all; 

% Init
temin=18.573; % keV
temax=18.58; % keV

nte=1000;
tex=5; % year
eres = 0.01e-3; % keV
mnu = 1e-3; % 1 eV
%mnu = 7; % keV

% KATRIN : tritium mass 
global com_opt1 ; com_opt1 = {...
    'tex',tex,'RNS_nTnu',1000,...
    'nTe',nte,'Temin',temin,...
    'Temax',temax,'Tmass',1000,'energy_resol',...
    eres,'mnu',mnu};
A1   = TritiumRelicNu(com_opt1{:},'Tmass',100); A1.RateCaptureT_KATRIN(1);
A10  = TritiumRelicNu(com_opt1{:},'Tmass',1000); A10.RateCaptureT_KATRIN(1);
A100 = TritiumRelicNu(com_opt1{:},'Tmass',100000); A100.RateCaptureT_KATRIN(1);

% Relic Neutrino Spectrum
figure(1)
h10 = semilogy((A1.Te-A1.Q)*1000,A10.T_eCapture_S,'LineWidth',2,'Color','Black','LineStyle','-');
hold on
%line([A1.mnu A1.mnu]*1000,[0 max(A1.T_eCapture_S)],'LineWidth',4,'Color','Red')
h1 = semilogy((A1.Te-A1.Q)*1000,A1.T_eCapture_S,'LineWidth',2,'Color','Blue','LineStyle','-.');
h100 = semilogy((A1.Te-A1.Q)*1000,A100.T_eCapture_S,'LineWidth',2,'Color','Red','LineStyle','--');
grid on;
xlabel('Kinetic Energy (keV)','FontSize',12);
str = sprintf('dN/dE (per %.1s keV)',A1.TeStep);
ylabel(str,'FontSize',14);
clear str1, clear str10, clear srt100;
str1 = sprintf('E-capture: Res.=%.1s eV , %.0f gT, %.0f evts',A1.sigma_E*1e3,100, A1.RateCaptureT_KATRIN(1));
str10 = sprintf('E-capture: Res.=%.1s eV , %.0f gT, %.0f evts',A10.sigma_E*1e3,1000,A10.RateCaptureT_KATRIN(1));
str100 = sprintf('E-capture: Res.=%.1s eV , %.0f gT, %.0f evts',A100.sigma_E*1e3,10000,A100.RateCaptureT_KATRIN(1));
PrettyFigureFormat;
%publish_figure(1,'NuCEspectrum.eps');

% Tritium Spectrum + Electron Capture
global com_opt2 ; com_opt2 = {...
    'NormType','RelicNu','Tmass',1,...
    'tex',tex,...
    'nTe',nte,'Temin',temin,'Temax',temax,...
    'mnu4',0,'sin22th4',0,...
    'mnu',1e-3,...
    'FermiFuncType',0,...
    'RadiativeC',0,...
    'RecoilCoulombC',0,...
    'RecoilWmVmAC',0,...
    'EEexchangeC',0,...
    'ScreeningC',0,...
    'WintFiniteSizeC',0,...
    'FiniteExtChargeC',0,...
    'T2exStates',0,...
    'energy_resol',1e-4
    };
clear T; 
T1=TritiumSpectrum(com_opt2{:},'Tmass',100);
T10=TritiumSpectrum(com_opt2{:},'Tmass',1000);
T100=TritiumSpectrum(com_opt2{:},'Tmass',10000);

hT1 = semilogy((T1.Te-A1.Q)*1000,T1.Spectrum,'LineWidth',2,'Color','Black','LineStyle','-');
hT10 = semilogy((T1.Te-A1.Q)*1000,T10.Spectrum,'LineWidth',2,'Color','Blue','LineStyle','-.');
hT100 = semilogy((T1.Te-A1.Q)*1000,T100.Spectrum,'LineWidth',2,'Color','Red','LineStyle','--');
grid on;
xlabel('E-E_0 (eV)','FontSize',12);
str = sprintf('dN/dE (per %.1s eV)',A1.TeStep*1000);
ylabel(str,'FontSize',14);
strT1 = sprintf('Tritium Beta Decay: Energy Res = %.1s eV',A1.sigma_E*1000);
lh1 = legend([h1 h10 h100 hT1],str1,str10,str100,strT1);
%legend(lh1,'box','off');
set(lh1,'FontSize',12);
%axis([-1 2 1e-2 10000])
title('Relic Neutrinos Capture Tritium - 5 years','FontSize',14);

PrettyFigureFormat;
hold off
publish_figure(1,'relic.eps');

return;

%% Tritium Spectrum With/Without Neutrino Mass
temin = 18.575-10e-3; 
temax = 18.575+1e-3;
global com_opt3 ; com_opt3 = {...
    'NormType','KATRIN',...
    'tex',tex,...
    'nTe',10000,...
    'mnu4',0,'sin22th4',0,...
    'FermiFuncType',0,...
    'RadiativeC',0,...
    'RecoilCoulombC',0,...
    'RecoilWmVmAC',0,...
    'EEexchangeC',0,...
    'ScreeningC',0,...
    'WintFiniteSizeC',0,...
    'FiniteExtChargeC',0,...
    'T2exStates',0,...
    'energy_resol',0,...
    };
Tnomass=TritiumSpectrum(com_opt3{:},'mnu',0);
T1eV=TritiumSpectrum(com_opt3{:},'mnu',1e-3);

hTnomass = semilogy((Tnomass.Te-Tnomass.Q)*1000,Tnomass.Spectrum,'LineWidth',2,'Color','Black','LineStyle','--');
hold on
hT1eV = semilogy((Tnomass.Te-Tnomass.Q)*1000,T1eV.Spectrum,'LineWidth',2,'Color','Blue','LineStyle','-');
hold off
grid on;
xlabel('E-Q (eV)','FontSize',12);
ylabel(str,'FontSize',14);
lh1 = legend([hTnomass hT1eV],'no mass','1 eV');
%legend(lh1,'boxoff');
set(lh1,'FontSize',12);
%axis([-30e-3 1e-3 1e-2 1e10])
title('KATRIN Tritium Spectrum - 5 years','FontSize',14);

PrettyFigureFormat;
hold off
publish_figure(1,'relic.eps');
