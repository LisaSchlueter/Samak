%% Plots for KSN1 sensitivity contours

%% Data Importation

% Each datafile contains three variables : 
%   sith4_X (mixing angle)
%   m4_Y    (sterile mass)
%   chi_Z   (chiSq values, if everything is OK this should be only 4.61)

filepath   = [getenv('SamakPath'),'ksn1ana/contour/contourmatfiles/'];
CLflag     = '95'; % Confidence Level
file_A     = sprintf('coord_95eV_Real_stat_%s_thierry.mat',CLflag);
file_B     = sprintf('coord_95eV_Real_syst_%s_thierry.mat',CLflag);

da  = importdata([filepath,file_A]);   % Data A
db  = importdata([filepath,file_B]);   % Data B

% Constant data
d_giunti   = importdata([filepath,'coord_Giunti.mat']);         % KATRIN Data from Giunti

d_raa_90   = importdata([filepath,'coord_RAA_90.mat']);         % Data RAA
d_raa_95_a = importdata([filepath,'coord_RAA_95_A.mat']);
d_raa_95_b = importdata([filepath,'coord_RAA_95_B.mat']);

d_mainz    = importdata([filepath,'coord_mainz.mat']);
d_troitsk  = importdata([filepath,'coord_troitsk.mat']);

d_stereo   = importdata([filepath,'coord_Stereo_95CL.mat']);

d_prospect = importdata([filepath,'coord_Prospect_95CL.mat']);

d_neutrino4 = importdata([filepath,'coord_Neutrino4_123sigma.mat']);

d_danss    = importdata([filepath,'coord_DANSS_95CL.mat']);

%% sin(2*th4)
da2  = 1-(1-2*da.sith4_X).^2;
db2  = 1-(1-2*db.sith4_X).^2;

%% Plot tunings
% Cutting the tails
na   = length(da2); nb   = length(db2);  
cuta = (3:na);  cutb = (3:nb);  

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
%title('KATRIN Sterile Neutrino Analysis (KS1) - 90% Sensitivity')
fig = figure('Renderer','painters');
set(fig, 'Units', 'normalized', 'Position', [0.001, 0.001,0.5, 0.999]);

% === Constant Data ===
%   Mainz & Troitsk
p_m     =       plot (smooth(d_mainz.sith4_X), smooth(d_mainz.m4_Y),...
                    '-.','color',rgb('Salmon'),'LineWidth',2);
hold on
p_t     =     plot (smooth(d_troitsk.sith4_X), smooth(d_troitsk.m4_Y),...
                    '-.','color',rgb('DarkSlateGrey'),'LineWidth',2);
hold on
%   RAA
p_raa95 =       plot (raa95x,raa95y,...
                    'color',rgb('ForestGreen')*1.2,'LineWidth',2);
hold on
p_raa2  =       plot (d_raa_95_b.sith4_X, d_raa_95_b.m4_Y,...
                    'color',rgb('ForestGreen')*1.2,'LineWidth',2);
hold on

%   Giunti
%p_g     =       plot (d_giunti.sith4_X, d_giunti.m4_Y,...
%                    'color',[0.9290 0.6940 0.1250],'LineWidth',1);
hold on
                
% Stereo 
p_s     =       plot (smooth(d_stereo.SinSquare2Theta_X),smooth(d_stereo.DmSquare41_Y),...
                    'color',rgb('Salmon'),'LineWidth',1);
                
% Prospect 
p_p     =       plot (smooth(d_prospect.SinSquare2Theta_X),smooth(d_prospect.DmSquare41_Y),...
                    'color',rgb('Violet'),'LineWidth',1);
                
% Neutrino4 
bn4     =       boundary((d_neutrino4.SinSquare2Theta_X_2sigma)',(d_neutrino4.DmSquare41_Y_2sigma)',0);
p_n     =       plot((d_neutrino4.SinSquare2Theta_X_2sigma(bn4)),(d_neutrino4.DmSquare41_Y_2sigma(bn4)),...
                    'color',rgb('Red'),'LineWidth',1);

% danss 
p_d     =       plot (smooth(d_danss.SinSquare2Theta_X),smooth(d_danss.DmSquare41_Y),...
                    'color',rgb('Grey'),'LineWidth',1);
                
% === Our Data ===
pA      =       plot (smooth([da2(cuta),1]),smooth([da.m4_Y(cuta),da.m4_Y(na)]),...
                    '--','color',prlB,'LineWidth',3);
hold on
pB      =       plot (smooth([db2(cutb),1]),smooth([db.m4_Y(cutb),db.m4_Y(nb)]),...
                    '-','color',prlB,'LineWidth',5);
%% Plot Parameters

% Axis
xlabel('sin^2(2\theta_{ee})');
ylabel('\Deltam_{41}^2  (eV^2)');

% Labels
mainz     = 'Mainz 90% CL';
troitsk   = 'Troitsk 90% CL';
stereo    = 'Stéréo 95% CL';
prospect  = 'Prospect 95% CL';
neutrino4 = 'Neutrino-4 2\sigma';
danss     = 'DANSS 95% CL';
katrinA   = ['KATRIN KSN1 - stat - ' CLflag '% CL'];% - [E_0-90;E_0+50] eV'];
katrinB   = ['KATRIN KSN1 - stat+syst - ' CLflag '% CL'];% - [E_0-90;E_0+50] eV'];
giunti    = 'arXiv:1912.12956 - - 90% CL - [E_0-40;E_0+50] eV';
raa90     = 'RAA+GA - PRD 83, 073006 (2011) - 90% CL';
raa95     = 'RAA+GA - PRD 83, 073006 (2011) - 95% CL';

hl=legend([p_m p_t p_s p_p p_n p_d pA pB p_raa95],...          % Label order omit p_g
        {mainz,troitsk,stereo,prospect,neutrino4,danss,...
        katrinA,katrinB,...                     % KATRIN
        raa95},...                              % RAA
        'Location','northoutside',...              % Legend settings
        'box','off')
%        giunti,...                              % Giunti
hl.NumColumns=2;
grid on

PRLFormat;

% Axis 
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');

axis([0.01 1 0.1 10000])
axis square
%title('KATRIN Sterile Neutrino Analysis (KSN1) - 95% CL Sensitivity') % Exclusion Limit (Data)
fileString = sprintf('./plots/ksn1_exlimit_%sCL_Dm24sinSq2T.pdf',CLflag);
export_fig(fig,fileString)