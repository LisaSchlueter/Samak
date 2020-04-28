%% Plots for KSN1 sensitivity contours with version 2+ files

%% Data Importation

% Each datafile contains these variables : 
% 'contour_settings','R','sith4_X','si2th4_X','m4_Y','m_beta','DM2','chi_Z','X0'

plt_title = 'KATRIN Sterile Neutrino Analysis (KS1) - 90% Sensitivity'

filepath   = [getenv('SamakPath'),'ksn1ana/contour/contour_files/V2'];
file_A     = 'coord_90eV_Real_syst_95.mat';
file_B     = 'coord_90eV_Real_syst_95_freeM.mat';
% file_C     = 'coord_90eV_Real_syst_99.mat';

da  = importdata([filepath,file_A]);    % Data A
setA = da.contour_settings;             % Contour settings A
db  = importdata([filepath,file_B]);    % Data B
setB = db.contour_settings;             % Contour settings B
% dc  = importdata([filepath,file_C]);  % Data c
% setC = dc.contour_settings;           % Contour settings C

% Constant data
d_giunti   = importdata([filepath,'coord_Giunti.mat']);         % KATRIN Data from Giunti

d_raa_90   = importdata([filepath,'coord_RAA_90.mat']);         % Data RAA
d_raa_95_a = importdata([filepath,'coord_RAA_95_A.mat']);
d_raa_95_b = importdata([filepath,'coord_RAA_95_B.mat']);

d_mainz    = importdata([filepath,'coord_mainz.mat']);          % Data Mainz
d_troitsk  = importdata([filepath,'coord_troitsk.mat']);        % Data Troitsk

%% Plot tunings
% Cutting the tails
na   = length(da.sith4_X); nb   = length(db.sith4_X); %nc   = length(dc.sith4_X);
cuta = (1:na);  cutb = (1:nb);  %cutc = (1:nc);

% Continuation of the RAA curves
raa90x = d_raa_90.sith4_X;   raa90y = d_raa_90.m4_Y;   n90 = length(raa90x);
raa95x = d_raa_95_a.sith4_X; raa95y = d_raa_95_a.m4_Y; n95 = length(raa95x);

raa90x = [raa90x(1),raa90x,raa90x(n90)];
raa95x = [raa95x(1),raa95x,raa95x(n95)];
raa95y = [10000,raa95y,10000]; raa90y = [10000,raa90y,10000];

%% Colors
prlG = [81 126 102]/255;
prlB = [50 148 216]/255;


%%%% ================ %%%%
%% Plotting
figure
title(plt_title)

% === Constant Data ===
%   Mainz & Troitsk
p_m     =       plot (d_mainz.sith4_X, d_mainz.m4_Y,...
                    'color',rgb('Gray'));
hold on
p_t     =       plot (d_troitsk.sith4_X, d_troitsk.m4_Y,...
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
pA      =       plot ([da.si2th4_X(cuta),1],[da.m4_Y(cuta),da.m4_Y(na)],...
                    'color',prlB,'LineWidth',3);
hold on

pB      =       plot ([db.si2th4_X(cutb),1],[db.m4_Y(cutb),db.m4_Y(nb)],...
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
giunti   = 'arXiv:1912.12956 - - 90%CL - [E_0-40;E_0+50] eV';
raa90    = 'Phys. Rev. D 83, 073006 (2011) - 90%CL';
raa95    = 'Phys. Rev. D 83, 073006 (2011) - 95%CL';

% KATRIN Labels

katrinA  = sprintf('KATRIN KSN1 %1$s %2$s - %3$s - %4$d % CL - [E_0-%5$d;E_0+50] eV',...
    s(setA.datatype),s(setA.activeFlag),s(setA.uncertainty),s(setA.CL),s(setA.eVrange));
katrinB  = sprintf('KATRIN KSN1 %1$s %2$s - %3$s - %4$d % CL - [E_0-%5$d;E_0+50] eV',...
    s(setB.datatype),s(setB.activeFlag),s(setB.uncertainty),s(setB.CL),s(setB.eVrange));
% katrinC  = sprintf('KATRIN KSN1 %1$s %2$s - %3$s - %4$d % CL - [E_0-%5$d;E_0+50] eV',...
%     s(setC.datatype),s(setC.activeFlag),s(setC.uncertainty),s(setC.CL),s(setC.eVrange));

% Legend
legend([p_m p_t pA pB p_g p_raa95 p_raa90],...      % Label order
        {mainz,troitsk,...                          % Mainz&Troitsk
        katrinA,katrinB,...                         % KATRIN
        giunti,...                                  % Giunti
        raa95,raa90},...                            % RAA
        'Location','southwest',...                  % Legend settings
        'box','off')

grid on

PRLFormat;

% Axis 
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');

axis([0.01 1 0.1 10000])
axis square

%% MISC
function b = s(a)
% Just a function for automatic labels
    if strcmp(a,'syst')
        b = 'stat+syst';
    elseif strcmp(a,'stat')
        b='stat';
    end
    if strcmp(a,'Twin')
        b='MCdata';
    elseif strcmp(a,'Real')
        b='data';
    end
    if strcmp(a,'FREE')
        b='free m_\beta';
    elseif strcmp(a,'FIX')
        b='m_\beta^2 = -1';
    elseif strcmp(a,'OFF')
        b='m_\beta = 0';
    end
end