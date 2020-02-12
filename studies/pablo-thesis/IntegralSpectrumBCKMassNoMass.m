clear; close all;

% Add path to s folder if it does not exist
addpath(genpath('../../../Samak2.0'));

% Spectrum with mnu = 0 eV (actually mnu^2 = 0 ev^2)
A0 = ref_sensitivity('nTeBinningFactor',100,'TD','Flat30');
A0.ComputeTBDDS();
A0.ComputeTBDIS();

% Spectrum with mnu = 1 eV (actually mnu^2 = 1 ev^2)
A1 = ref_sensitivity('nTeBinningFactor',100,'TD','Flat30');
A1.mnuSq_i = (0.2)^2;
A1.ComputeTBDDS();
A1.ComputeTBDIS();

% Spectrum with mnu = 0 eV (actually mnu^2 = 0 ev^2)
B0 = ref_sensitivity('nTeBinningFactor',100,'TD','Flat30','BKG_RateAllFPDSec',0.4);
B0.ComputeTBDDS();
B0.ComputeTBDIS();

% Spectrum with mnu = 1 eV (actually mnu^2 = 1 ev^2)
B1 = ref_sensitivity('nTeBinningFactor',100,'TD','Flat30','BKG_RateAllFPDSec',0.4);
B1.mnuSq_i = (0.2)^2;
B1.ComputeTBDDS();
B1.ComputeTBDIS();

lw = 2;

SnuS0 = A1.TBDIS./A0.TBDIS;
SnuS0B = B1.TBDIS./B0.TBDIS;

hold on
h2 = bar(A1.qU - A1.Q_i, SnuS0-1,'FaceColor',rgb('LightGray'),'EdgeColor',rgb('DarkGray'));
% h3 = errorbar(A1.qU - A1.Q_i, SnuS0-1, A1.TBDISE./A0.TBDIS,'.','Color',rgb('DarkSlateGrey'),'LineWidth',1);

h4 = bar(A1.qU - A1.Q_i, SnuS0B-1,'FaceColor',rgb('Tomato'),'EdgeColor',rgb('IndianRed'));
% h5 = errorbar(A1.qU - A1.Q_i, SnuS0B-1, B1.TBDISE./B0.TBDIS,'.','Color',rgb('DarkRed'),'LineWidth',1);
hold off

title('Tritium \beta-decay spectra')
ylabel('S_{\nu}/S_0-1');
xlabel('E - E_0 (eV)');

legend([h2,h4],{'BCK = 10 mcps','BCK = 400 mcps'},'Location','Best');

% limits on axis to show endpoint
%axis([-3,+0.5,0,0.001])


PrettyFigureFormat;
