%Comparison to Formaggio: 
%effective Tmass = 66.5e-6g
%Measurement time 1 yr
%expected events measured: 5e-6

ToggleES = 'ON';
mnuSq_i = 0;
eta_i = 1;
RhoD = (66.5e-6*0.57)/(pi*4.5^2*1e3*0.95*0.083*(2*3.0160492e-3/6.022140857e23));


% init TBD object with Design Report (TDR) properties:
T1 = ref_RelicNuBkg_DesignReport('ToggleES',ToggleES,'mnuSq_i',mnuSq_i,'eta_i',eta_i,'TimeSec',365.242*24*3600);
T2 = ref_RelicNuBkg_DesignReport('ToggleRelic','OFF','mnuSq_i',mnuSq_i,'TimeSec',365.242*24*3600);

% Compute diff spec.
T1.ComputeTBDDS;
T2.ComputeTBDDS;
T1.TBDDS = T1.TBDDS*RhoD/T1.WGTS_CD_MolPerCm2;
T2.TBDDS = T2.TBDDS*RhoD/T2.WGTS_CD_MolPerCm2;

fig1=figure(1);
hT1 = semilogy((T1.Te-T1.Q),T1.TBDDS,'LineWidth',2,'Color','Black','LineStyle','-');
hold on;
hT2 = semilogy((T1.Te-T1.Q),T2.TBDDS,'LineWidth',2,'Color','Blue','LineStyle',':');
hA1 = semilogy((T1.Te-T1.Q),T1.TBDDS_R,'LineWidth',2,'Color','Red','LineStyle','--');
grid on;
ylim([1e-20,1e2]);
xlabel('E-E_0 (eV)','FontSize',12);
str = sprintf('dN/dE');
ylabel(str,'FontSize',14);
strT2 = sprintf('Tritium Beta Decay Spectrum');
strT1 = sprintf('Beta Decay + E-capture');
str1 = sprintf('E-capture: %.1g evts',T1.NormFactorTBDDS_R*RhoD/T1.WGTS_CD_MolPerCm2);
lh1 = legend([hT1 hT2 hA1],strT1,strT2,str1);
%legend(lh1,'box','off');
set(lh1,'FontSize',12);
%axis([-1 2 1e-2 10000])
title(sprintf('Relic Neutrinos Capture Tritium - %.1f years',T1.TimeSec/(365.242*24*3600)),'FontSize',14);
PrettyFigureFormat;
hold off;

% ini TBD object with Design Report (TDR) properties:
A = ref_RelicNuBkg_DesignReport('DopplerEffectFlag','FSD');
A.DE_sigma = 0.0935;
A.LoadFSD;
% Compute diff spec.
A.ComputeTBDDS;
A.TBDDS = A.TBDDS./simpsons(A.Te,A.TBDDS);
pFSD = A.PlotTBDDS;
hold on;

A.TTFSD = 'OFF';
A.DTFSD = 'OFF';
A.HTFSD = 'OFF';

A.ComputeTBDDS;
pNoFSD = plot(A.Te-A.Q,A.TBDDS./simpsons(A.Te,A.TBDDS),'-.','LineWidth',3,'Color',rgb('Orange'));
hold off;
%%
leg = legend([pNoFSD,pFSD],'Without FSD','With FSD');
leg.Title.String = sprintf('Differential \\beta-spectrum');
leg.Title.FontWeight = 'normal';
xlim([-15 1]);
set(gca,'YScale','log');
ylim([5e-100 10]);
