close all; clear all;
addpath('../../tools/masumhabib-PlotPub-23bcfed/lib/');
mtd      = '6hMTD';
runtimes = 86400;
nfit     = 1000;

WGTS_TASR_RelErr_V = 0.01*[0 1 2 3 4 5 6 7 8 9 10]';

% Activity Fluctuation
muhat_tasr    = [-0.3462      -1.3154      -0.7976      -1.3425     -10.8871       3.6372      -6.0636       3.1383     -12.7686     -10.6034      -0.8415];
sigmahat_tasr = [180.7411      223.4369      276.7608       315.813      366.2993      393.3558      437.4058      471.5154      503.6088      549.7951      600.3166];

% Response Function
muhat_rf    = ones(1,11).*(-5.6);
sigmahat_rf = ones(1,11).*(337.3);


% Activiy Flucuation + Response Function
muhat_tasrrf    = [-5.6         19.3        -25.7         59.4         55.8        -20.4           37        -22.1         53.2         19.8         17.9];
sigmahat_tasrrf = [337.3          402.7          465.9          563.5          686.7          764.5          823.3          882.7          962.3           1027         1119.5];


figure
plot(WGTS_TASR_RelErr_V*100,sigmahat_tasr);
opt = [];
opt.LineStyle = {'-.'}; 
opt.Markers = {'s'}; 
opt.MarkerSpacing = [ 1];
pltstr     = sprintf('(E0)eff - 6h MTD - 1 day measurement');
opt.Title  = pltstr; % plot title
opt.XLabel = 'Subrun Activity Fluctuation within a run (%)'; % xlabel
opt.YLabel = '(E0)_{eff} uncertainty (meV)'; %ylabel
opt.YScale = 'lin'; % 'linear' or 'log'
opt.XScale = 'lin'; % 'linear' or 'log'
opt.FontSize = 16;
opt.XLim       = [0 4];
opt.YLim       = [0 800];
opt.Legend = {'stat. + T-activity fluctuation'};
opt.LegendLoc = 'NorthWest';
plttitle   = sprintf('./results/E0eff_%s_%gs_%gfits_TASR.png',mtd,runtimes,nfit);
opt.FileName = plttitle; 
setPlotProp(opt);


figure
plot(WGTS_TASR_RelErr_V*100,sigmahat_rf);
opt = [];
opt.LineStyle = {'--' }; 
opt.Colors = [1,0,0]; 
opt.Markers = { 'd' }; 
opt.MarkerSpacing = [1 ];
pltstr     = sprintf('(E0)eff - 6h MTD - 1 day measurement');
opt.Title  = pltstr; % plot title
opt.XLabel = 'Subrun Activity Fluctuation within a run (%)'; % xlabel
opt.YLabel = '(E0)_{eff} uncertainty (meV)'; %ylabel
opt.YScale = 'lin'; % 'linear' or 'log'
opt.XScale = 'lin'; % 'linear' or 'log'
opt.FontSize = 16;
opt.XLim       = [0 4];
opt.YLim       = [0 800];
opt.Legend = {'stat + RF'};
opt.LegendLoc = 'NorthWest';
plttitle   = sprintf('./results/E0eff_%s_%gs_%gfits_RF.png',mtd,runtimes,nfit);
opt.FileName = plttitle; 
setPlotProp(opt);

figure
plot(WGTS_TASR_RelErr_V*100,sigmahat_tasr,WGTS_TASR_RelErr_V*100,sigmahat_rf);
opt = [];
opt.LineStyle = {'-.','--'}; 
opt.Markers = {'s','d'}; 
opt.MarkerSpacing = [1 1 1];
pltstr     = sprintf('(E0)eff - 6h MTD - 1 day measurement');
opt.Title  = pltstr; % plot title
opt.XLabel = 'Subrun Activity Fluctuation within a run (%)'; % xlabel
opt.YLabel = '(E0)_{eff} uncertainty (meV)'; %ylabel
opt.YScale = 'lin'; % 'linear' or 'log'
opt.XScale = 'lin'; % 'linear' or 'log'
opt.FontSize = 16;
opt.XLim       = [0 4];
opt.YLim       = [0 800];
opt.Legend = {'stat. + T-activity fluctuation','stat + RF'};
opt.LegendLoc = 'NorthWest';
plttitle   = sprintf('./results/E0eff_%s_%gs_%gfits_TASRandRF.png',mtd,runtimes,nfit);
opt.FileName = plttitle; 
setPlotProp(opt);

figure
plot(WGTS_TASR_RelErr_V*100,sigmahat_tasr,WGTS_TASR_RelErr_V*100,sigmahat_rf,WGTS_TASR_RelErr_V*100,sigmahat_tasrrf);
opt = [];
opt.LineStyle = {'-.','--','-'}; 
opt.Markers = {'s','d','o'}; 
opt.MarkerSpacing = [1 1 1];
pltstr     = sprintf('(E0)eff - 6h MTD - 1 day measurement');
opt.Title  = pltstr; % plot title
opt.XLabel = 'Subrun Activity Fluctuation within a run (%)'; % xlabel
opt.YLabel = '(E0)_{eff} uncertainty (meV)'; %ylabel
opt.YScale = 'lin'; % 'linear' or 'log'
opt.XScale = 'lin'; % 'linear' or 'log'
opt.FontSize = 16;
opt.XLim       = [0 4];
opt.YLim       = [0 800];
opt.Legend = {'stat. + T-activity fluctuation','stat + RF','stat. + T-activity fluctuation + RF'};
opt.LegendLoc = 'NorthWest';
plttitle   = sprintf('./results/E0eff_%s_%gs_%gfits_TASR+RF.png',mtd,runtimes,nfit);
opt.FileName = plttitle; 
setPlotProp(opt);
