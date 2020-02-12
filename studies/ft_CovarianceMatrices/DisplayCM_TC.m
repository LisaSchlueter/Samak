% Init
TD = 'Flat60';
TimeSec =900*24*60*60;
MACE_Ba_T = 7*1e-04;
MACE_Bmax_T = 6*0.7;
WGTS_B_T    = 3.6*0.7;
WGTS_CD_MolPerCm2 = 5e17;
SysEffect = 'TC';
nTrials = 1000;

A= ref_TBD_NominalKATRIN('TD',TD,'TimeSec',TimeSec,'recomputeRF','OFF',...
    'MACE_Ba_T',MACE_Ba_T,'WGTS_B_T',WGTS_B_T,'MACE_Bmax_T',MACE_Bmax_T,...
    'WGTS_CD_MolPerCm2',WGTS_CD_MolPerCm2);
A.FiniteExtChargeFlag = 'OFF';
A.WintFiniteSizeFlag  = 'OFF';
A.EEexchangeFlag      = 'OFF';
A.RecoilCoulombFlag   = 'OFF';
A.RecoilWmVmAFlag     = 'OFF';
A.RadiativeFlag       = 'OFF';
A.ComputeTBDDS; A.ComputeTBDIS;
TBDIS_NC = A.TBDIS; % Without any Corrections
TBDDS_NC = A.TBDDS;
A.FiniteExtChargeFlag = 'ON';
A.WintFiniteSizeFlag  = 'ON';
A.EEexchangeFlag      = 'ON';
A.RecoilCoulombFlag   = 'ON';
A.RecoilWmVmAFlag     = 'ON';
A.RadiativeFlag       = 'ON';
A.ComputeTBDDS; A.ComputeTBDIS;
TBDIS_TC = A.TBDIS; % With all Corrections
TBDDS_TC = A.TBDDS;

%%
A.FiniteExtChargeFlag = 'ON';
A.WintFiniteSizeFlag  = 'OFF';
A.EEexchangeFlag      = 'OFF';
A.RecoilCoulombFlag   = 'OFF';
A.RecoilWmVmAFlag     = 'OFF';
A.RadiativeFlag       = 'OFF';
A.ComputeTBDDS; A.ComputeTBDIS;
TBDIS_TC1 = A.TBDIS; 
TBDDS_TC1 = A.TBDDS;

A.FiniteExtChargeFlag = 'OFF';
A.WintFiniteSizeFlag  = 'ON';
A.EEexchangeFlag      = 'OFF';
A.RecoilCoulombFlag   = 'OFF';
A.RecoilWmVmAFlag     = 'OFF';
A.RadiativeFlag       = 'OFF';
A.ComputeTBDDS; A.ComputeTBDIS;
TBDIS_TC2 = A.TBDIS; 
TBDDS_TC2 = A.TBDDS;

A.FiniteExtChargeFlag = 'OFF';
A.WintFiniteSizeFlag  = 'OFF';
A.EEexchangeFlag      = 'ON';
A.RecoilCoulombFlag   = 'OFF';
A.RecoilWmVmAFlag     = 'OFF';
A.RadiativeFlag       = 'OFF';
A.ComputeTBDDS; A.ComputeTBDIS;
TBDIS_TC3 = A.TBDIS; 
TBDDS_TC3 = A.TBDDS;

A.FiniteExtChargeFlag = 'OFF';
A.WintFiniteSizeFlag  = 'OFF';
A.EEexchangeFlag      = 'OFF';
A.RecoilCoulombFlag   = 'ON';
A.RecoilWmVmAFlag     = 'OFF';
A.RadiativeFlag       = 'OFF';
A.ComputeTBDDS; A.ComputeTBDIS;
TBDIS_TC4 = A.TBDIS; 
TBDDS_TC4 = A.TBDDS;

A.FiniteExtChargeFlag = 'OFF';
A.WintFiniteSizeFlag  = 'OFF';
A.EEexchangeFlag      = 'OFF';
A.RecoilCoulombFlag   = 'OFF';
A.RecoilWmVmAFlag     = 'ON';
A.RadiativeFlag       = 'OFF';
A.ComputeTBDDS; A.ComputeTBDIS;
TBDIS_TC5 = A.TBDIS; 
TBDDS_TC5 = A.TBDDS;

A.FiniteExtChargeFlag = 'OFF';
A.WintFiniteSizeFlag  = 'OFF';
A.EEexchangeFlag      = 'OFF';
A.RecoilCoulombFlag   = 'OFF';
A.RecoilWmVmAFlag     = 'OFF';
A.RadiativeFlag       = 'ON';
A.ComputeTBDDS; A.ComputeTBDIS;
TBDIS_TC6 = A.TBDIS; 
TBDDS_TC6 = A.TBDDS;
%% plot integral spectrum
f3 = figure('Renderer','opengl');
set(f3, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7 ,0.5]);
qUlong = linspace(A.qU(1),A.qU(end),1000);
plot(qUlong-18575,interp1(A.qU,TBDIS_TC,qUlong,'spline')./interp1(A.qU,TBDIS_NC,qUlong,'spline'),'Color',rgb('CadetBlue'),'LineWidth',4);
xlabel('retarding potential qU - 18575 (V)');
ylabel('S^{TC} / S^{NC}');
PrettyFigureFormat;
set(gca,'FontSize',20);
grid on;
xlim([A.qU(1)-18575 A.qU(end)-18575]);
ylim([0.999995*min(TBDIS_TC./TBDIS_NC) 1.00001]);
f3_str = sprintf('TC_IS-Ratio');
publish_figurePDF(f3,['./plots/CovMatInfo/pdf/',f3_str,'.pdf']);
print(f3,['./plots/CovMatInfo/png/',f3_str,'.png'],'-dpng');
savefig(f3,['./plots/CovMatInfo/fig/',f3_str,'.fig'],'compact');

%% Plot differential spectrum
close
f3 = figure('Renderer','opengl');
set(f3, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7 ,0.7]);
plot(A.Te-18575,TBDDS_TC./TBDDS_NC,'Color',rgb('CadetBlue'),'LineWidth',4);
%plot(A.qU,TBDIS_TC./TBDIS_NC,'Color',rgb('CadetBlue'),'LineWidth',4);
xlabel('kin. energy - 18575 (eV)');
ylabel('\Gamma^{TC} / \Gamma^{NC}');
PrettyFigureFormat;
set(gca,'FontSize',20);

f3_str = sprintf('TC_DS-Ratio');
publish_figurePDF(f3,['./plots/CovMatInfo/pdf/',f3_str,'.pdf']);
print(f3,['./plots/CovMatInfo/png/',f3_str,'.png'],'-dpng');
savefig(f3,['./plots/CovMatInfo/fig/',f3_str,'.fig'],'compact');


%% plot all integral spectra
close
f3 = figure('Renderer','opengl');
set(f3, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7 ,0.7]);
qUlong = linspace(A.qU(1),A.qU(end),1000);
p1   = plot(qUlong-18575,interp1(A.qU,TBDIS_TC1,qUlong,'spline')./interp1(A.qU,TBDIS_NC,qUlong,'spline'),'Color',rgb('CadetBlue'),'LineWidth',2);
hold on;
p2   = plot(qUlong-18575,interp1(A.qU,TBDIS_TC2,qUlong,'spline')./interp1(A.qU,TBDIS_NC,qUlong,'spline'),'Color',rgb('RoyalBlue'),'LineWidth',2);
p3   = plot(qUlong-18575,interp1(A.qU,TBDIS_TC3,qUlong,'spline')./interp1(A.qU,TBDIS_NC,qUlong,'spline'),'Color',rgb('IndianRed'),'LineWidth',2);
p4   = plot(qUlong-18575,interp1(A.qU,TBDIS_TC4,qUlong,'spline')./interp1(A.qU,TBDIS_NC,qUlong,'spline'),'Color',rgb('GoldenRod'),'LineWidth',2);
p5   = plot(qUlong-18575,interp1(A.qU,TBDIS_TC5,qUlong,'spline')./interp1(A.qU,TBDIS_NC,qUlong,'spline'),'Color',rgb('DarkRed'),'LineWidth',2);
%p6   = plot(qUlong-18575,interp1(A.qU,TBDIS_TC6,qUlong,'spline')./interp1(A.qU,TBDIS_NC,qUlong,'spline'),'Color',rgb('Orange'),'LineWidth',2);
%pall = plot(qUlong-18575,interp1(A.qU,TBDIS_TC,qUlong,'spline')./interp1(A.qU,TBDIS_NC,qUlong,'spline'),'Color',rgb('IndianRed'),'LineWidth',2);
xlabel('retarding potential qU - 18575 (eV)');
ylabel('S^{TC} / S^{NC}');
PrettyFigureFormat;
set(gca,'FontSize',20);
grid on;
xlim([A.qU(1)-18575 A.qU(end)-18575]);
%ylim([0.999995*min(TBDIS_TC./TBDIS_NC) 1.00001]);
%cellstr(num2str(round(log10(get(gca,'YTick')), '10^%d'));
leg_str = {'finite extension of nucleus charge','weak interaction finite size correction','electron-electron exchange term',...
'recoil coulomb effect','^3He recoil corrections'};   
leg = legend(leg_str{:});
leg.FontSize = 16; leg.Location = 'best';
legend boxoff
f3_str = sprintf('TC_IS-Ratio');
publish_figurePDF(f3,['../ft_CovarianceMatrices/plots/CovMatInfo/pdf/',f3_str,'.pdf']);
print(f3,['../ft_CovarianceMatrices/plots/CovMatInfo/png/',f3_str,'.png'],'-dpng');
savefig(f3,['../ft_CovarianceMatrices/plots/CovMatInfo/fig/',f3_str,'.fig'],'compact');
