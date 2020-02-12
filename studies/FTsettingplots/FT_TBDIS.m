% First Tritium Data
% Plot TBDIS
%

close all

% Inputs
MTD=load('FT-TL2.mat');
TD      = MTD.TD;
TimeSec = MTD.RunTime*1;
CD      = 2.5e17;

%% All Pixels, 1 run
% All Pixels - Transmission function 50%
s50 = ref_FTMTD('WGTS_CD_MolPerCm2',CD,'TimeSec',TimeSec,'TD',TD);
s50.ComputeTBDDS;s50.ComputeTBDIS;

figure(1);
errorbar(s50.qU-s50.Q,s50.TBDIS,s50.TBDISE);
plt = Plot();
plt.LineWidth = 2; 
plt.LineStyle = '--'; 
plt.Markers = {'s'}; 
pltstr     = sprintf('Integral Spectrum All Pixels - %g sec - %s',TimeSec,TD);
plt.Title  = pltstr; % plot title
plt.XLabel = 'E-qU (V)'; % xlabel
plt.YLabel = 'Counts'; %ylabel
plt.YScale = 'log'; % 'linear' or 'log'
%plt.YLim       = [0 1.1];
plt.XScale = 'lin'; % 'linear' or 'log'
plt.FontSize = 16;
plt.Legend = {'\rho d = 2.5e17 / cm^2'}
plttitle   = sprintf('FT-TBDIS-AllPixel-%g-%s.png',TimeSec,TD);
plt.export(plttitle);


%% Single Pixel - Transmission function 50%
s50 = ref_FTMTD('WGTS_CD_MolPerCm2',CD,'TimeSec',TimeSec/148,'TD',TD);
s50.ComputeTBDDS;s50.ComputeTBDIS;

figure(2);
errorbar(s50.qU-s50.Q,s50.TBDIS,s50.TBDISE);
plt = Plot();
plt.LineWidth = 2; 
plt.LineStyle = '--'; 
plt.Markers = {'s'}; 
pltstr     = sprintf('Integral Spectrum 1 pixel - %g sec - %s',TimeSec,TD);
plt.Title  = pltstr; % plot title
plt.XLabel = 'E-qU (V)'; % xlabel
plt.YLabel = 'Counts'; %ylabel
plt.YScale = 'log'; % 'linear' or 'log'
%plt.YLim       = [0 1.1];
plt.XScale = 'lin'; % 'linear' or 'log'
plt.FontSize = 16;
plt.Legend = {'\rho d = 2.5e17 / cm^2'}
plttitle   = sprintf('FT-TBDIS-OnePixel-%g-%s.png',TimeSec,TD);
plt.export(plttitle);
