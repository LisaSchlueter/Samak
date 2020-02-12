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
Bkg       = [10 50 100 200 300 400 500 1000]; % mcps
stat30    = ([0.1800    0.2070    0.2220    0.2500    0.2690    0.2860    0.2990    0.3500] );
statsys30 =  ([0.2070    0.2260    0.2380    0.2620    0.2780    0.2940    0.3060    0.3540] );
stat60    =  ([0.1480    0.1770    0.1950    0.2120    0.2260    0.2410    0.2520    0.2930] );
statsys60 =  ([0.1880    0.2050    0.2170    0.2300    0.2410    0.2540    0.2630    0.3000] );

% sys
figure(1)
plot(Bkg(1:7),stat30(1:7),...
    Bkg(1:7),statsys30(1:7),...
    Bkg(1:7),stat60(1:7),...
    Bkg(1:7),statsys60(1:7));
plt = Plot();
plt.Colors = { % three colors for three data set
    [0.25, 0.25, 0.25] % data set 1
    [0.25, 0.25, 0.25] % data set 2
    [1, 0, 0] % data set 3
    [1, 0, 0]};
plt.YGrid = 'on'; % 'on' or 'off'
plt.BoxDim = [14, 6]; %[width, height] in inches
plt.Legend = {'30 eV Optimized MTD - stat', '30 eV Optimized MTD - stat+sys','60 eV Optimized MTD - stat', '60 eV Optimized MTD - stat+sys'}; % legends
plt.LegendLoc  =  'SouthEast';
plt.XLabel = 'Background (mcps)'; % xlabel
plt.YLabel = 'm_\nu Upper Limit - 90% CL'; %ylabel
plt.Title  = 'KATRIN stat + sys (0.017 eV^2) - 900 days - B_A=3 G Bs=3.6 T Bm=9 T - Samak Simulation'; % plot title
plt.LineWidth = [2, 5, 2, 5]; % three line widths
plt.LineStyle = {'--', '-','--', '-'}; % three line styles
plt.MarkerSpacing = [10, 100];
plt.YScale = 'lin'; % 'linear' or 'log'
plt.XScale = 'lin'; % 'linear' or 'log'
plt.FontSize = 24;
plt.export('Plot_KATRIN_900d_TDRsyst.png');

%% Background versus Ba
figure(2)
Ba             = [2:0.5:25]; % cps
BkgBa          = GetBackground('MACE_Ba_T',Ba*1e-4,'WGTS_B_T',3.6)*1e3;
BkgBa70        = GetBackground('MACE_Ba_T',Ba*1e-4,'WGTS_B_T',3.6*0.7)*1e3;
plot(Ba,BkgBa,Ba,BkgBa70);

[X,Y] = meshgrid(Ba,BkgBa);
Z = 18575*X./(6)*1e-4;
Z70 = 18575*X./(6*0.7)*1e-4;
hold on
[C,h]=contour(X,Y,Z,'ShowText','on');
hcl=clabel(C,h,'FontSize',15,'Color','red');
[C70,h70]=contour(X,Y,Z70,'ShowText','on');
clabel(C70,h70,'FontSize',15,'Color','blue');
hold off

plt = Plot();
plt.Colors = { % three colors for three data set
    rgb('red') , ...
    rgb('blue') , ...
    rgb('red') , ...
    rgb('blue') };
plt.YGrid = 'on'; % 'on' or 'off'
plt.BoxDim = [14, 6]; %[width, height] in inches
plt.Legend = {'B_m = 6 T','Bm = 4.5 T'}; % legends
plt.LegendLoc  =  'NorthEast';
plt.XLabel = 'Ba (G)'; % xlabel
plt.YLabel = 'Background (mcps)'; %ylabel
plt.Title  = 'KATRIN '; % plot title
plt.LineWidth = [5 5 2 2]; % three line widths
plt.LineStyle = {'--','-','-','--'}; % three line styles
plt.MarkerSpacing = [10, 100, 10 , 10];
plt.YScale = 'lin'; % 'linear' or 'log'
plt.XScale = 'lin'; % 'linear' or 'log'
plt.XLim   = [1 15];
plt.FontSize = 24;
plt.LegendBox = 'on';
plt.LegendBoxColor = rgb('white');
plt.export('Plot_KATRIN_Ba_Bkg.png');

%% Resolution (eV) versus Ba
figure(3)
Bs            = 6;
Bs70          = 4.5;
Ba            = [2:0.1:15]; % G
R             = 18575*Ba/Bs*1e-4;
R70           = 18575*Ba/Bs70*1e-4;
plot(Ba,R,Ba,R70);
plt = Plot();
plt.Colors = { % three colors for three data set
    [0.25, 0.25, 0.25] % data set 1
    [0.25, 0.25, 0.25]};
plt.YGrid = 'on'; % 'on' or 'off'
plt.BoxDim = [14, 6]; %[width, height] in inches
plt.Legend = {'Bm = 6 T','Bm = 4.5 T (70% nominal)'}; % legends
plt.LegendLoc  =  'SouthEast';
plt.XLabel = 'Ba (G)'; % xlabel
plt.YLabel = 'Resolution (eV)'; %ylabel
plt.Title  = 'KATRIN - MacE Filter Energy Resolution at E=18575 eV'; % plot title
plt.LineWidth = [4 4]; % three line widths
plt.LineStyle = {'--','-'}; % three line styles
plt.MarkerSpacing = [10, 10];
plt.YScale = 'lin'; % 'linear' or 'log'
plt.XScale = 'lin'; % 'linear' or 'log'
plt.FontSize = 30;
plt.export('Plot_KATRIN_Ba_Resolution.png');

%% Resolution (eV) versus Background
figure(4)
Bs            = 6;
Bs70          = 4.5;
Ba            = [2:0.1:15]; % G
R             = 18575*Ba/Bs*1e-4;
R70           = 18575*Ba/Bs70*1e-4;
BkgBa         = GetBackground('MACE_Ba_T',Ba*1e-4,'WGTS_B_T',3.6*0.7)*1e3;
plot(R70,BkgBa);
plt = Plot();
plt.Colors = { % three colors for three data set
    [0.25, 0.25, 0.25]};
plt.YGrid = 'on'; % 'on' or 'off'
plt.BoxDim = [14, 6]; %[width, height] in inches
plt.Legend = {'R = 18575 * Ba/Bm eV'}; % legends
plt.LegendLoc  =  'NorthEast';
plt.YLabel = 'Background (mcps)'; % xlabel
plt.XLabel = 'Resolution (eV)'; %ylabel
plt.Title  = 'KATRIN - MacE Resolution at E=18575 eV Versus Background'; % plot title
plt.LineWidth = [4]; % three line widths
plt.LineStyle = {'-'}; % three line styles
plt.MarkerSpacing = [10];
plt.YScale = 'lin'; % 'linear' or 'log'
plt.XScale = 'lin'; % 'linear' or 'log'
plt.FontSize = 30;
plt.XLim = [0.9 5];
plt.export('Plot_KATRIN_Background_Resolution.png');

%%

figure(1111)
Ba            = 3:1:12;          % G
Bs            = (3.6:0.1:6)/6*3.6;           % T 
Bm            = Bs*6/3.6;
[X,Y] = meshgrid(Ba,Bm);
range=60;
xxxSense = zeros(numel(Ba),numel(Bs));
i=0;
for ba=Ba
    j=0;i=i+1;
    for bs=Bs
        j=j+1;
        %TD = sprintf('Sensitivity_%.0feV_Ba%.0fG',range,ba);
        commonOpt = {'range',range,'mnuSq_i_Fit',(0:21:4000)*1e-03,'plotFit','ON','SysEffect','all',};
        [~,~,~,~,~,l90,mNumin] = ...
            NuMassScan_SensitivityNominal(...
            'TimeSec',179.5*86400*5,...
            commonOpt{:},'chi2','chi2Stat',...
            'WGTS_B_T',bs,'MACE_Ba_T',ba*1e-04,'MACE_Bmax_T',6*bs/3.6,...
            'plotFit','ON');
        xxxSense(i,j) = l90;
    end
end

R             = 18575*X./Y*1e-4; %eV
%%
B=GetBackground('MACE_Ba_T',X*1e-4,'WGTS_B_T',3.6.*(Y./6))*1e3; % mcps
figure(888)
[C,h]=contourf(X*1e-4,Y,sqrt(xxxSense)'); 
hcl=clabel(C,h,'FontSize',12,'Color','blue');
%colorbar; 
colormap(flipud(pink));
   hold on
   %[C1,h1]=contour(X*1e-4,Y,B,'ShowText','on');
   % hcl=clabel(C1,h1,'FontSize',20,'Color','black');
   %  [C2,h2]=contour(X*1e-4,Y,R,'ShowText','on');
   %  hcl=clabel(C2,h2,'FontSize',20,'Color','red');
   hold off

plt = Plot();
plt.YGrid = 'on'; % 'on' or 'off'
plt.XGrid = 'on'; % 'on' or 'off'
plt.BoxDim = [8, 8]; %[width, height] in inches
plt.XLabel = 'Ba (G)'; % xlabel
plt.YLabel = 'Bm (T)'; %ylabel
plt.Title  = 'KATRIN 5 x 180 days - 30eV MTD - stat'; % plot title
%plt.LineWidth = [4 4]; % three line widths
%plt.LineStyle = {'--','-'}; % three line styles
%plt.MarkerSpacing = [10, 10];
plt.YScale = 'lin'; % 'linear' or 'log'
plt.XScale = 'lin'; % 'linear' or 'log'
plt.FontSize = 30;
plt.export('xxx.png');