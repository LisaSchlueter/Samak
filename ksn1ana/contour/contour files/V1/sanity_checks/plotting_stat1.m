%% Plots for KSN1 sensitivity contours

%% Data Importation

% Each datafile contains three variables : 
%   sith4_X (mixing angle)
%   m4_Y    (sterile mass)
%   chi_Z   (chiSq values, if everything is OK this should be only 4.61)

filepath   = [getenv('SamakPath'),'ksn1ana/contour/sanity_checks/'];
file_A     = 'coord_stat+FPDeff.mat';
file_B     = 'coord_stat+FSD.mat';
file_C     = 'coord_stat+RF_BF.mat';
file_D     = 'coord_stat+RF_BX.mat';
file_E     = 'coord_stat+RF_EL.mat';
file_F     = 'coord_stat+Stack.mat';
file_G     = 'coord_stat+TASR.mat';
file_H     = 'coord_stat+TCoff_OTHER.mat';
file_I     = 'coord_stat+BKG.mat';

da  = importdata([filepath,file_A]);   % Data A
db  = importdata([filepath,file_B]);   % Data B
dc  = importdata([filepath,file_C]);   % Data C
dd  = importdata([filepath,file_D]);   % Data D
de  = importdata([filepath,file_E]);   % Data E
df  = importdata([filepath,file_F]);   % Data F
dg  = importdata([filepath,file_G]);   % Data G
dh  = importdata([filepath,file_H]);   % Data H
di  = importdata([filepath,file_I]);   % Data I

% Constant data
filepath   = [getenv('SamakPath'),'ksn1ana/contour/'];
d_giunti   = importdata([filepath,'coord_Giunti.mat']);         % KATRIN Data from Giunti

d_raa_90   = importdata([filepath,'coord_RAA_90.mat']);         % Data RAA
d_raa_95_a = importdata([filepath,'coord_RAA_95_A.mat']);
d_raa_95_b = importdata([filepath,'coord_RAA_95_B.mat']);

%% sin(2*th4)

da2  = 1-(1-2*da.sith4_X).^2;
db2  = 1-(1-2*db.sith4_X).^2;
dc2  = 1-(1-2*dc.sith4_X).^2;
dd2  = 1-(1-2*dd.sith4_X).^2;
de2  = 1-(1-2*de.sith4_X).^2;
df2  = 1-(1-2*df.sith4_X).^2;
dg2  = 1-(1-2*dg.sith4_X).^2;
dh2  = 1-(1-2*dh.sith4_X).^2;
di2  = 1-(1-2*di.sith4_X).^2;

%% Plot tunings

% Continuation of the RAA curves
raa90x = d_raa_90.sith4_X;   raa90y = d_raa_90.m4_Y;   n90 = length(raa90x);
raa95x = d_raa_95_a.sith4_X; raa95y = d_raa_95_a.m4_Y; n95 = length(raa95x);

raa90x = [raa90x(1),raa90x,raa90x(n90)];
raa95x = [raa95x(1),raa95x,raa95x(n95)];
raa95y = [10000,raa95y,10000]; raa90y = [10000,raa90y,10000];

%% Colors
prlG = [81 126 102]/255;
prlB = [50 148 216]/255;

%% Plotting
figure

% === Constant Data ===
%   RAA
p_raa90 =       plot (raa90x,raa90y,...
                    'color',prlG,'LineWidth',1);
hold on
p_raa95 =       plot (raa95x,raa95y,...
                    'color',rgb('ForestGreen')*1.2,'LineWidth',1);
hold on
p_raa2  =       plot (d_raa_95_b.sith4_X, d_raa_95_b.m4_Y,...
                    'color',rgb('ForestGreen')*1.2,'LineWidth',1);
hold on

%   Giunti
p_g     =       plot (d_giunti.sith4_X, d_giunti.m4_Y,...
                    'color',[0.9290 0.6940 0.1250],'LineWidth',1);

                
% === Our Data ===
strt = 6;
pA      =       plot (da2(strt:end),da.m4_Y(strt:end),...
                    'color',rgb('yellow'),'LineWidth',3);
hold on
pB      =       plot (db2(strt:end),db.m4_Y(strt:end),...
                    '--','color',rgb('red'),'LineWidth',3);
hold on
pC      =       plot (dc2(strt:end),dc.m4_Y(strt:end),...
                    '--','color',rgb('blue'),'LineWidth',3);
hold on
pD      =       plot (dd2(strt:end),dd.m4_Y(strt:end),...
                    '--','color',rgb('green'),'LineWidth',3);
hold on
pE      =       plot (de2(strt:end),de.m4_Y(strt:end),...
                    ':','color',rgb('orange'),'LineWidth',3);
hold on
pF      =       plot (df2(strt:end),df.m4_Y(strt:end),...
                    '--','color',rgb('black'),'LineWidth',3);
hold on
pG      =       plot (dg2(strt:end),dg.m4_Y(strt:end),...
                    '--','color',rgb('DeepPink'),'LineWidth',3);
hold on
pH      =       plot (dh2(strt:end),dh.m4_Y(strt:end),...
                    '--','color',rgb('Cyan'),'LineWidth',3);
hold on
pI      =       plot (di2(strt:end),di.m4_Y(strt:end),...
                    '--','color',rgb('DarkSlateGray'),'LineWidth',3);


%% Plot Parameters

% Axis
xlabel('sin^2(2\theta_{ee})');
ylabel('\Deltam_{41}^2  (eV^2)');

% Labels
katrinA  = 'FPD';
katrinB  = 'Final State Distribution';
katrinC  = 'Magnetic Fields';
katrinD  = 'Column Density';
katrinE  = 'Energy Loss function';
katrinF  = 'HV Fluctuations';
katrinG  = 'Tritium activity fluctuations';
katrinH  = 'Theoretical corrections';
katrinI  = 'Background';

giunti   = 'arXiv:1912.12956 - - 90%CL - [E_0-40;E_0+50] eV';
raa90    = 'Phys. Rev. D 83, 073006 (2011) - 90%CL';
raa95    = 'Phys. Rev. D 83, 073006 (2011) - 95%CL';

legend([pA pB pC pD pE pF pG pH pI p_g p_raa95 p_raa90],...     % Label order
        {katrinA,katrinB,katrinC,katrinD,katrinE,...
        katrinF,katrinG,katrinH,katrinI...                      % KATRIN
        giunti,...                                              % Giunti
        raa95,raa90},...                                        % RAA
        'Location','southwest',...                              % Legend settings
        'box','off')

grid on

PRLFormat;

% Axis 
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');

axis([0.01 1 0.1 10000])
axis square

title('KATRIN Sterile Neutrino Analysis (KSN1) - 90% CL Exclusion Limit, statistics only +1 systematic') % Exclusion Limit (Data)