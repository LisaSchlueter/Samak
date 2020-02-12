% Ba
Ba = 3:.1:15; MACE_Ba_T = Ba*1e-4;

% FT
MACE_Ba_T_RelErr_FT = 2.*MACE_Ba_T./MACE_Ba_T;

% CM35
MACE_Ba_T_RelErrCM35 = 1.6*1e-06./MACE_Ba_T*100; %fixed to absolute value of 2*10^-6 T

% SS
MACE_Ba_T_RelErr_SS = (MACE_Ba_T*3.9*1e-03+1.07*1e-06)./MACE_Ba_T*100;  % Nikolaus


% Plot
figure(1)
plt = Plot(Ba,MACE_Ba_T_RelErr_FT,Ba,MACE_Ba_T_RelErr_SS,Ba,MACE_Ba_T_RelErrCM35);
plt.Colors = { rgb('DarkBlue') ,rgb('CadetBlue') , rgb('FireBrick')};
plt.YGrid = 'on';  
plt.YLim = [1e-1,5]; 
plt.BoxDim = [14, 6]; %[width, height] in inches
plt.Legend = {'2% - used for FT','N. Trost - used for sensitivity studies','constant? - KATRIN CM35'}; % legends
plt.LegendLoc  =  'NorthEast';
plt.XLabel = 'Ba (G)'; % xlabel
plt.YLabel = 'Ba Relative Uncertainty (%)'; %ylabel
plt.Title  = 'AP Magnetic Field Uncertainty (Ba)'; % plot title
plt.LineWidth = [8, 8 , 8]; % three line widths
plt.LineStyle = {'-', '-' , '-'}; % three line styles
plt.MarkerSpacing = [10, 100];
plt.YScale = 'log'; % 'linear' or 'log'
plt.XScale = 'lin'; % 'linear' or 'log'
plt.FontSize = 24;
plt.export('BaUncertainties.png');

%%
figure(2)
plt = Plot(Ba,MACE_Ba_T_RelErr_FT.*Ba/100*18575/0.7/3.6e4*1e3,Ba,MACE_Ba_T_RelErr_SS.*Ba/100*18575/0.7/3.6e4*1e3,Ba,MACE_Ba_T_RelErrCM35.*Ba/100*18575/0.7/3.6e4*1e3);
plt.Colors = { rgb('DarkBlue') ,rgb('CadetBlue') , rgb('FireBrick')};
plt.YGrid = 'on';  
plt.YLim = [1e-2,0.3]*1e3; 
plt.BoxDim = [14, 6]; %[width, height] in inches
plt.Legend = {'2% - used for FT','N. Trost - used for sensitivity studies','constant? - KATRIN CM35'}; % legends
plt.LegendLoc  =  'NorthWest';
plt.XLabel = 'Ba (G)'; % xlabel
plt.YLabel = 'Uncertainty (meV) '; %ylabel
plt.Title  = 'AP B-field induced E-resolution Uncertainty'; % plot title
plt.LineWidth = [8, 8 , 8]; % three line widths
plt.LineStyle = {'-', '-' , '-'}; % three line styles
plt.MarkerSpacing = [10, 100];
plt.YScale = 'log'; % 'linear' or 'log'
plt.XScale = 'lin'; % 'linear' or 'log'
plt.FontSize = 24;
plt.export('BaUncertaintiesDeltaReseV.png');
