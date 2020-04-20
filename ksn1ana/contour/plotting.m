%% Plots for KSN1 sensitivity contours

%% Data Importation

% Each datafile contains three variables : 
%   sith4_X (mixing angle)
%   m4_Y    (sterile mass)
%   chi_Z   (chiSq values, if everything is OK this should be only 4.61)

filepath   = [getenv('SamakPath'),'ksn1ana/contour/'];
file_A     = 'coord_fitrium_40_Twin_95.mat';
file_B     = 'coord_90eV_Twin_stat_95_newN.mat';
% file_C     = 'coord_90eV_Real_syst_99.mat';

da  = importdata([filepath,file_A]);   % Data A
db  = importdata([filepath,file_B]);   % Data B
% dc  = importdata([filepath,file_C]);   % Data c

% Constant data
d_giunti   = importdata([filepath,'coord_Giunti.mat']);         % KATRIN Data from Giunti

d_raa_90   = importdata([filepath,'coord_RAA_90.mat']);         % Data RAA
d_raa_95_a = importdata([filepath,'coord_RAA_95_A.mat']);
d_raa_95_b = importdata([filepath,'coord_RAA_95_B.mat']);

d_mainz    = importdata([filepath,'coord_mainz.mat']);
d_troitsk  = importdata([filepath,'coord_troitsk.mat']);

%% sin(2*th4)

da2  = 1-(1-2*da.sith4_X).^2;
db2  = 1-(1-2*db.sith4_X).^2;
% dc2  = 1-(1-2*dc.sith4_X).^2;

%% Plot tunings
% Cutting the tails
na   = length(da2); nb   = length(db2); %nc   = length(dc2);
cuta = (13:na);  cutb = (4:nb);  %cutc = (5:nc);

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
title('KATRIN Sterile Neutrino Analysis (KS1) - 90% Sensitivity')

% === Constant Data ===
%   Mainz & Troitsk
p_m     =       plot (d_mainz.sith4_X, d_mainz.m4_Y,...
                    'color',rgb('Gray'));
hold on
p_t     =     plot (d_troitsk.sith4_X, d_troitsk.m4_Y,...
                    'color',rgb('Navy'));
hold on
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
hold on
                
% === Our Data ===
pA      =       plot ([da2(cuta),1],[da.m4_Y(cuta),da.m4_Y(na)],...
                    'color',prlB,'LineWidth',3);
hold on
pB      =       plot ([db2(cutb),1],[db.m4_Y(cutb),db.m4_Y(nb)],...
                    '--','color',prlB,'LineWidth',3);
% hold on
% pC      =       plot ([dc2(cutc),1],[dc.m4_Y(cutc),dc.m4_Y(nc)],...
%                     ':','color',prlB,'LineWidth',3);

%% Plot Parameters

% Axis
xlabel('sin^2(2\theta_{ee})');
ylabel('\Deltam_{41}^2  (eV^2)');

% Labels
mainz    = 'Mainz data 90%CL';
troitsk  = 'Troitsk data 90%CL';
katrinA  = 'KATRIN fitrium MCdata - stat - 95%CL - [E_0-90;E_0+50] eV';
katrinB  = 'KATRIN KSN1 MCdata - stat - 90%CL - [E_0-90;E_0+50] eV';
katrinC  = 'KATRIN KSN1 data - stat+sys - 99%CL - [E_0-90;E_0+50] eV';
giunti   = 'arXiv:1912.12956 - - 90%CL - [E_0-40;E_0+50] eV';
raa90    = 'Phys. Rev. D 83, 073006 (2011) - 90%CL';
raa95    = 'Phys. Rev. D 83, 073006 (2011) - 95%CL';

legend([p_m p_t pA pB p_g p_raa95 p_raa90],...          % Label order
        {mainz,troitsk,...
        katrinA,katrinB,...                    % KATRIN
        giunti,...                              % Giunti
        raa95,raa90},...                        % RAA
        'Location','southwest',...              % Legend settings
        'box','off')

grid on

PRLFormat;

% Axis 
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');

axis([0.01 1 0.1 10000])
axis square

title('KATRIN Sterile Neutrino Analysis (KSN1) - 95% CL Sensitivity') % Exclusion Limit (Data)