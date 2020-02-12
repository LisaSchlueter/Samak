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

lw = 2;

SnuS0 = A1.TBDIS./A0.TBDIS;

hold on
h2 = bar(A1.qU - A1.Q_i, SnuS0-1,'FaceColor',rgb('LightGray'),'EdgeColor',rgb('DarkGray'));
h3 = errorbar(A1.qU - A1.Q_i, SnuS0-1, A1.TBDISE./A0.TBDIS,'.','Color',rgb('DarkSlateGrey'),'LineWidth',1);
% bar(A1.qU-A1.Q_i, SnuS0-1,'LineWidth',lw);
hold off

title('Tritium \beta-decay spectra [mass / no mass]')
ylabel('S_{\nu}/S_0-1');
xlabel('E - E_0 (eV)');

% limits on axis to show endpoint
%axis([-3,+0.5,0,0.001])


% arrows with annotations 
% fs = 16;
% x0 = [0.6 0.41];
% y0 = [0.7 0.41];
% a0 = annotation('textarrow',x0,y0,'String','m_{\nu} = 0 eV');
% a0.FontSize = fs;
% 
% x1 = [0.31 0.35];
% y1 = [0.25 0.35];
% a1 = annotation('textarrow',x1,y1,'String','m_{\nu} = 1 eV');  
% a1.FontSize = fs;

PrettyFigureFormat;

% export_fig('plots/SpecWithMass.pdf','-pdf')
% copyfile plots/SpecWithMass.pdf ../../../../Nextcloud/master-thesis/LatexThesis/figs/
