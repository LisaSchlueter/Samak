ToggleES = 'ON';
mnuSq_i = 1;
eta_i = 1e11;

CD = 1./(2*pi*4.5^2*0.973); %column density: one tritium atom
CD = 5e17;

% init TBD object with Design Report (TDR) properties:
T1 = ref_RelicNuBkg_TDR('ToggleES',ToggleES,'mnuSq_i',mnuSq_i,'eta_i',eta_i,'WGTS_CD_MolPerCm2',CD);
T2 = ref_RelicNuBkg_TDR('ToggleRelic','OFF','mnuSq_i',mnuSq_i,'WGTS_CD_MolPerCm2',CD);

%compute diff spec.
T1.ComputeTBDDS;
T2.ComputeTBDDS;
%T1.PlotTBDDS;
%T2.PlotTBDDS;

% Compute int spec.
T1.ComputeTBDIS;
T2.ComputeTBDIS;
T1.TBDIS = T1.TBDIS./T1.qUfrac;
T2.TBDIS = T2.TBDIS./T2.qUfrac;
T1.PlotTBDIS('type','log');
T2.PlotTBDIS('type','log');

% hT1 = semilogy((T1.Te-T1.Q),T1.TBDIS,'LineWidth',2,'Color','Black','LineStyle','-');
% hold on;
% hT2 = semilogy((T1.Te-T1.Q),T2.TBDIS,'LineWidth',2,'Color','Blue','LineStyle',':');
% hA1 = semilogy((T1.Te-T1.Q),T1.TBDIS_R,'LineWidth',2,'Color','Red','LineStyle','--');
% grid on;
% ylim([1e-20,1e2]);
% xlabel('E-E_0 (eV)','FontSize',12);
% str = sprintf('dN/dE');
% ylabel(str,'FontSize',14);
% strT1 = sprintf('Tritium Beta Decay Spectrum');
% strT2 = sprintf('Beta Decay + E-capture');
% str1 = sprintf('E-capture: %.1g evts',T1.NormFactorTBDIS_R);
% lh1 = legend([hT1 hT2 hA1],strT1,strT2,str1);
% %legend(lh1,'box','off');
% set(lh1,'FontSize',12);
% %axis([-1 2 1e-2 10000])
% title(sprintf('Relic Neutrinos Capture Tritium - %.1f years',T1.TimeSec/(365.242*24*3600)),'FontSize',14);
% PrettyFigureFormat;
% hold off;