% First Tritium
% Plot Response Function expected
% Thierry, 2/5/18

close all
% RF Display
qUtest = 18500;
EqUmin = 0;
EqUmax = 1999;

% Transmission function 25%
s25 = ref_FTMTD('WGTS_CD_MolPerCm2',1.25e17);
s25.ComputeTBDDS;s25.ComputeTBDIS;
% Transmission function 50%
s50 = ref_FTMTD('WGTS_CD_MolPerCm2',2.50e17);
s50.ComputeTBDDS;s50.ComputeTBDIS;
% Transmission function 100%
s100 = ref_FTMTD('WGTS_CD_MolPerCm2',5.00e17);
s100.ComputeTBDDS;s100.ComputeTBDIS;

%% Plot 1
e  = qUtest+EqUmin:0.01:qUtest+EqUmax;
plt = Plot(e-qUtest,s25.KTF(e,qUtest),e-qUtest,s50.KTF(e,qUtest),e-qUtest,s100.KTF(e,qUtest));
plt.Title = 'Response Function - qU=18500 V'; % plot title
plt.XLabel = 'E-qU (V)'; % xlabel
plt.YLabel = 'Transmission Probability'; %ylabel
plt.YScale = 'lin'; % 'linear' or 'log'
plt.YLim       = [0 1.1];
plt.XScale = 'log'; % 'linear' or 'log'
plt.FontSize = 16;
plt.Legend = {'\rho d = 1.25e17 cm^{-2}', '\rho d = 2.5e17 cm^{-2}', '\rho d = 5e17 cm^{-2}'}; % legends
plt.LegendLoc = 'SouthEast';
plt.export('FT-RFs.png');

%% Plot 2: Zoom
e  = qUtest-1:0.01:qUtest+5;
plt = Plot(e-qUtest,s25.KTF(e,qUtest),e-qUtest,s50.KTF(e,qUtest),e-qUtest,s100.KTF(e,qUtest));
plt.Title = 'Response Function - qU=18500 V'; % plot title
plt.XLabel = 'E-qU (V)'; % xlabel
plt.YLabel = 'Transmission Probability'; %ylabel
plt.YScale = 'lin'; % 'linear' or 'log'
plt.YLim       = [-0.1 0.9];
plt.XLim       = [-0.5 4];
plt.XScale = 'lin'; % 'linear' or 'log'
plt.FontSize = 16;
plt.Legend = {'\rho d = 1.25e17 cm^{-2}', '\rho d = 2.5e17 cm^{-2}', '\rho d = 5e17 cm^{-2}'}; % legends
plt.LegendLoc = 'NorthWest';
plt.export('FT-RFs-zoom.png');

%% Plot 3
clear plt;
plt = Plot(s25.qU-s25.Q,s25.TBDIS,s25.qU-s25.Q,s50.TBDIS,s25.qU-s25.Q,s100.TBDIS);
plt.Title = 'TBDIS - Whole FPD'; % plot title
plt.XLabel = 'E-qU (V)'; % xlabel
pltstr = sprintf('counts for a %g sec run',s25.TimeSec);
plt.YLabel = pltstr; %ylabel
plt.YScale = 'log'; % 'linear' or 'log'
%plt.YLim       = [0 1.1];
plt.XScale = 'lin'; % 'linear' or 'log'
plt.FontSize = 16;
plt.Legend = {'\rho d = 1.25e17 cm^{-2}', '\rho d = 2.5e17 cm^{-2}', '\rho d = 5e17 cm^{-2}'}; % legends
plt.LegendLoc = 'NorthEast';
plt.export('FT-TDBISs.png');

%% Plot 3
clear plt;
plt = Plot(s25.qU-s25.Q,s25.TBDIS./s100.TBDIS,s25.qU-s25.Q,s50.TBDIS./s100.TBDIS,s25.qU-s25.Q,s100.TBDIS./s100.TBDIS);
plt.Title = 'TBDIS'; % plot title
plt.XLabel = 'E-qU (V)'; % xlabel
plt.YLabel = 'ratio'; %ylabel
plt.YScale = 'lin'; % 'linear' or 'log'
%plt.YLim       = [0 1.1];
plt.XScale = 'lin'; % 'linear' or 'log'
plt.FontSize = 16;
plt.Legend = {'\rho d = 1.25e17 cm^{-2}', '\rho d = 2.5e17 cm^{-2}', '\rho d = 5e17 cm^{-2}'}; % legends
plt.LegendLoc = 'NorthWest';
plt.export('FT-TDBISs-ratio.png');