%
% Plot KATRIN TDR Sensitivity 
% as a function of backgrounds
% 
% Ba=3G
% Bs=3.6T
% Bm=6T
%
% stat: 900 days
% sys:0.017 eV^2
%
% Optimized MTD 3G
% - 30 eV
% - 60 eV
% 

% data
Bkg         =  [50 100 200 364 ]; % mcps
mtd30       =  ( [254 270 294 319] );
mtd30_nosfd =  ( [254 270 294 319] );
mtd60       =  ( [242 260 283 306] );
mtd60_nosfd =  ( [229 244 264 286] );
av          = (mtd30 + mtd60_nosfd)/2;

% sys
figure(1)
plot(Bkg,mtd30,...
    Bkg,mtd60_nosfd,...
    Bkg,av);
plt = Plot();
plt.Colors = { % three colors for three data set
    rgb('DarkGreen') % data set 1
    rgb('IndianRed')
    rgb('DarkBlue')};
plt.YGrid = 'on'; % 'on' or 'off'
plt.BoxDim = [14, 6]; %[width, height] in inches
plt.Legend = {'30 eV MTD - Convervative','60 eV MTD - No FSD & No Eloss systematics - Challenging','Recommandation for Design Studies'}; % legends
plt.LegendLoc  =  'SouthEast';
plt.XLabel = 'Background (mcps)'; % xlabel
plt.YLabel = 'm_\nu Upper Limit - 90% CL'; %ylabel
plt.Title  = 'KATRIN 1260 days (stat + sys) - B_A=7 G Bs=2.5 T Bm=4.2 T (Samak)'; % plot title
plt.LineWidth = [7,7, 4]; % three line widths
plt.LineStyle = {'-', '-', '--'}; % three line styles
plt.MarkerSpacing = [10, 10, 10];
plt.YScale = 'lin'; % 'linear' or 'log'
plt.XScale = 'lin'; % 'linear' or 'log'
plt.FontSize = 24;
plt.AxisLineWidth = 4;
plt.TickLength = [0.01, 0.01];
plt.export('Plot_KATRIN_1260d_Background.png');