%% Plots for KSN1 sensitivity contours

%% Data Importation

% Each datafile contains three variables : 
%   sith4_X (mixing angle)
%   m4_Y    (sterile mass)
%   chi_Z   (chiSq values, if everything is OK this should be only 4.61)

filepath       = [getenv('SamakPath'),'ksn1ana/contour/sanity_checks/'];
filepathstat   = [getenv('SamakPath'),'ksn1ana/contour/contourmatfiles/'];
file_A     = 'coord_stat+FPDeff.mat';
file_B     = 'coord_stat+FSD.mat';
file_C     = 'coord_stat+RF_BF.mat';
file_D     = 'coord_stat+RF_BX.mat';
file_E     = 'coord_stat+RF_EL.mat';
file_F     = 'coord_stat+Stack.mat';
file_G     = 'coord_stat+TASR.mat';
file_H     = 'coord_stat+TCoff_OTHER.mat';
file_I     = 'coord_stat+BKG.mat';
file_J     = 'coord_90eV_Real_stat.mat';

da  = importdata([filepath,file_A]);   % Data A
db  = importdata([filepath,file_B]);   % Data B
dc  = importdata([filepath,file_C]);   % Data C
dd  = importdata([filepath,file_D]);   % Data D
de  = importdata([filepath,file_E]);   % Data E
df  = importdata([filepath,file_F]);   % Data F
dg  = importdata([filepath,file_G]);   % Data G
dh  = importdata([filepath,file_H]);   % Data H
di  = importdata([filepath,file_I]);   % Data I
dj  = importdata([filepathstat,file_J]);   % Data J


%% sin(2*th4)

da2  = da.sith4_X;
db2  = db.sith4_X;
dc2  = dc.sith4_X;
dd2  = dd.sith4_X;
de2  = de.sith4_X;
df2  = df.sith4_X;
dg2  = dg.sith4_X;
dh2  = dh.sith4_X;
di2  = di.sith4_X;
dj2  = dj.sith4_X;

%% Colors
prlG = [81 126 102]/255;
prlB = [50 148 216]/255;

%% Plotting
fig = figure('Renderer','painters');
set(fig, 'Units', 'normalized', 'Position', [0.001, 0.001,0.5, 0.999]);
                
% === Our Data ===
strt = 6;
pA      =       plot ((da2(strt:end)),(sqrt(da.m4_Y(strt:end))),...
                    'color',rgb('yellow'),'LineWidth',4);
hold on
pB      =       plot ((db2(strt:end)),(sqrt(db.m4_Y(strt:end))),...
                    '--','color',rgb('red'),'LineWidth',4);
hold on
pC      =       plot ((dc2(strt:end)),(sqrt(dc.m4_Y(strt:end))),...
                    '--','color',rgb('blue'),'LineWidth',4);
hold on
pD      =       plot ((dd2(strt:end)),(sqrt(dd.m4_Y(strt:end))),...
                    '--','color',rgb('green'),'LineWidth',4);
hold on
pE      =       plot ((de2(strt:end)),(sqrt(de.m4_Y(strt:end))),...
                    ':','color',rgb('orange'),'LineWidth',4);
hold on
pF      =       plot ((df2(strt:end)),(sqrt(df.m4_Y(strt:end))),...
                    '--','color',rgb('black'),'LineWidth',4);
hold on
pG      =       plot ((dg2(strt:end)),(sqrt(dg.m4_Y(strt:end))),...
                    '--','color',rgb('DeepPink'),'LineWidth',4);
hold on
pH      =       plot ((dh2(strt:end)),(sqrt(dh.m4_Y(strt:end))),...
                    '--','color',rgb('Cyan'),'LineWidth',4);
hold on
pI      =       plot ((di2(strt:end)),(sqrt(di.m4_Y(strt:end))),...
                    '--','color',rgb('DarkSlateGray'),'LineWidth',4);
hold on
pJ      =       plot ((dj2(strt:end)),(sqrt(dj.m4_Y(strt:end))),...
                    '-','color',rgb('Black'),'LineWidth',4);
                

%% Plot Parameters

% Axis
xlabel('|U_{e4}|^2');
ylabel('m_{4}  (eV)');

% Labels
katrinA  = 'Stat + FPD';
katrinB  = 'Stat + Final State Distribution';
katrinC  = 'Stat + Magnetic Fields';
katrinD  = 'Stat + Column Density';
katrinE  = 'Stat + Energy Loss function';
katrinF  = 'Stat + HV Fluctuations';
katrinG  = 'Stat + Tritium activity fluctuations';
katrinH  = 'Stat + Theoretical corrections';
katrinI  = 'Stat + Background Slope';
katrinJ  = 'KSN1 Stat. Only - 90%CL';

legend([pA pB pC pD pE pF pG pH pI pJ],...     % Label order
        {katrinA,katrinB,katrinC,katrinD,katrinE,...
        katrinF,katrinG,katrinH,katrinI,katrinJ...                      % KATRIN
        },...                                        % RAA
        'Location','southwest',...                              % Legend settings
        'box','off')

grid on

PRLFormat;

% Axis 
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');


axis([0.005 0.5 0.5 100])
axis square
export_fig(fig,'./plots/ksn1_contour_test3.pdf')

%title('KATRIN Sterile Neutrino Analysis (KSN1) - 90% CL Exclusion Limit, statistics only +1 systematic') % Exclusion Limit (Data)