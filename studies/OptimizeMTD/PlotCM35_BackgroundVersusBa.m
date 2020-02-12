%% Background versus Ba
figure(2)
Ba               = [2:0.2:12]; % cps
BkgBa            = GetBackground('MACE_Ba_T',Ba*1e-4,'WGTS_B_T',3.6,'Anchor6G','OFF')*1e3;
BkgBa70          = GetBackground('MACE_Ba_T',Ba*1e-4,'WGTS_B_T',3.6*0.7,'Anchor6G','OFF')*1e3;
BkgBa70FTroi26        = GetBackground('MACE_Ba_T',Ba*1e-4,'WGTS_B_T',3.6*0.7,'Anchor6G','ON','Anchor6GValue',335e-3)*1e3;
BkgBa70FTroi14        = GetBackground('MACE_Ba_T',Ba*1e-4,'WGTS_B_T',3.6*0.7,'Anchor6G','ON','Anchor6GValue',430e-3)*1e3;
plot(Ba,BkgBa70,Ba,BkgBa70FTroi26,Ba,BkgBa70FTroi14);

[X,Y] = meshgrid(Ba,BkgBa70);
%Z = 18575*X./(6)*1e-4;
Z70 = 18575*X./(6*0.7)*1e-4;
hold on
%[C,h]=contour(X,Y,Z,'ShowText','on');
%hcl=clabel(C,h,'FontSize',15,'Color','red');
[C70,h70]=contour(X,Y,Z70,'ShowText','on');
clabel(C70,h70,'FontSize',20,'Color',rgb('FireBrick'));
hold off

plt = Plot();
plt.Colors = { % three colors for three data set
    rgb('CadetBlue') , ...
    rgb('DarkBlue') , ...
    rgb('DarkBlue') , ...
    rgb('DarkBlue') , ...
    rgb('Black'), ...
    rgb('Black')};
plt.YGrid = 'on'; % 'on' or 'off'
plt.BoxDim = [14, 6]; %[width, height] in inches
plt.Legend = {'Neutrino 2016 scaled for B_m = 4.2 T , B_s = 2.5 T','First Tritium ROI [26-32] keV - B_m = 4.2 T , B_s = 2.5 T','First Tritium ROI [14-32] keV - B_m = 4.2 T , B_s = 2.5 T','MAC-E Resolution (eV)''MAC-E Resolution (eV)'}; % legends
plt.LegendLoc  =  'NorthEast';
plt.XLabel = 'Ba (G)'; % xlabel
plt.YLabel = 'Background (mcps)'; %ylabel
%plt.Title  = 'KATRIN Background Verus Anaysis Plane B-field '; % plot title
plt.LineWidth = [5 8 2 1]; % three line widths
plt.LineStyle = {'--','-','-.','--','--'}; % three line styles
plt.MarkerSpacing = [10, 100, 10 , 10, 10];
plt.YScale = 'lin'; % 'linear' or 'log'
plt.XScale = 'lin'; % 'linear' or 'log'
plt.XLim   = [2 11];
plt.YLim   = [100 800];
plt.FontSize = 24;
plt.LegendBox = 'on';
plt.LegendBoxColor = rgb('white');
%plt.export('Plot_KATRIN_Ba_Bkg.png');


