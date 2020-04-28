%% Plots for KSN1 sensitivity contours
%% In m4 sin2theta Plane

RealTwinFlag = 'Real';
MainzTroitzk = 'ON';

%% Data Importation

% Each datafile contains three variables : 
%   sith4_X (mixing angle)
%   m4_Y    (sterile mass)
%   chi_Z   (chiSq values, if everything is OK this should be only 4.61)

filepath   = [getenv('SamakPath'),'ksn1ana/contour/contourmatfiles/'];
CLflag     = '95'; % Confidence Level

switch RealTwinFlag
    case 'Real'
        TwinLabel   = '';
        % zero fixed nu mass - stat
        file_A1     = sprintf('coord_95eV_Real_stat_%s_thierry.mat',CLflag);
        % zero fixed nu mass - stat+sys
        file_A2     = sprintf('coord_90eV_Real_syst_%s.mat',CLflag);
%        file_A2     = sprintf('coord_95eV_Real_syst_%s_thierry.mat',CLflag);
        % free nu mass - stat+sys
        file_B      = sprintf('coord_90eV_Real_syst_%s_freeM.mat',CLflag);
%        file_B      = sprintf('coord_95eV_Real_syst_%s_thierry.mat',CLflag);
        % minus one fixed nu mass - stat+sys
        file_C      = sprintf('coord_90eV_Real_syst_%s.mat',CLflag);
%        file_C      = sprintf('coord_95eV_Real_syst_%s_thierry.mat',CLflag);
    case 'Twin'
        TwinLabel   = 'Simulation';
        % zero fixed nu mass - stat
        file_A1     = 'coord_90eV_Twin_stat_95.mat';
        % zero fixed nu mass - stat+sys
        file_A2     = 'coord_90eV_Twin_syst.mat';
        % free nu mass - stat+sys
        file_B      = 'coord_90eV_Twin_stat_95_freeM.mat';
        % minus one fixed nu mass - stat+sys
        file_C      = 'coord_90eV_Twin_stat_95_freeM.mat';
end

da1  = importdata([filepath,file_A1]);   % Data A1
da2  = importdata([filepath,file_A2]);   % Data A2
db   = importdata([filepath,file_B]);   % Data B
dc   = importdata([filepath,file_C]);   % Data C

% Constant data
d_giunti   = importdata([filepath,'coord_Giunti.mat']);         % KATRIN Data from Giunti
d_mainz    = importdata([filepath,'coord_mainz.mat']);
d_troitsk  = importdata([filepath,'coord_troitsk.mat']);

%% sin(th4)
da12  = da1.sith4_X;
da22  = da2.sith4_X;
db2   = db.sith4_X;
dc2   = dc.sith4_X;

%% Plot tunings
% Cutting the tails
na1   = length(da12); na2   = length(da22); nb   = length(db2);  nc   = length(dc2);  
cuta1 = (1:na1); cuta2 = (1:na2);  cutb = (1:nb);  cutc = (1:nc);  

%% Colors
prlA = [50 148 216]/255;
prlB = rgb('IndianRed');
prlC = rgb('DarkGreen');

%% Plotting
%title('KATRIN Sterile Neutrino Analysis (KS1) - 90% Sensitivity')
fig = figure('Renderer','painters');
set(fig, 'Units', 'normalized', 'Position', [0.001, 0.001,0.5, 0.999]);

%hold on
%   Giunti
%p_g     =       plot (d_giunti.sith4_X, d_giunti.m4_Y,...
%                    'color',[0.9290 0.6940 0.1250],'LineWidth',1);
                
hold on

% === Our Data ===
pA1      =       plot (smooth([da12(cuta1),1]),smooth(sqrt([da1.m4_Y(cuta1),da1.m4_Y(na1)])),...
    '--','color',prlA,'LineWidth',4);

pA2      =       plot (smooth([da22(cuta2),1]),smooth(sqrt([da2.m4_Y(cuta2),da2.m4_Y(na2)])),...
    '-','color',prlA,'LineWidth',4);

hold on
pB      =       plot (smooth([db2(cutb),1]),smooth(sqrt([db.m4_Y(cutb),db.m4_Y(nb)])),...
    '--','color',prlB,'LineWidth',4);

switch RealTwinFlag
    case 'Real'
        hold on
        pC      =       plot (smooth([dc2(cutc),1]),smooth(sqrt([dc.m4_Y(cutc),dc.m4_Y(nc)])),...
            '-.','color',prlC,'LineWidth',4);
end

switch MainzTroitzk
    case 'ON'
        
        %   Mainz & Troitsk
        mainzSinTheta2 = 0.5 * (1-sqrt(1-d_mainz.sith4_X));
        p_m     =       plot (smooth(mainzSinTheta2), smooth(sqrt(d_mainz.m4_Y)),...
            '-.','color',rgb('Salmon'),'LineWidth',2);
        % hold on
        troitskSinTheta2 = 0.5 * (1-sqrt(1-d_troitsk.sith4_X));
        p_t     =      plot (smooth(troitskSinTheta2), smooth(sqrt(d_troitsk.m4_Y)),...
            '-.','color',rgb('DarkSlateGrey'),'LineWidth',2);
end

%% Plot Parameters

% Axis
xlabel('|U_{e4}|^2');
ylabel('m_{4}  (eV)');

% Labels
mainz    = 'Mainz 90%CL - m_\nu^2 = 0 eV^2';
troitsk  = 'Troitsk 90%CL - m_\nu^2 = 0 eV^2';
katrinA1  = ['KATRIN KSN1 ' CLflag '%CL' TwinLabel ' - stat only - m_\nu^2 = 0 eV^2'];
katrinA2  = ['KATRIN KSN1 ' CLflag '%CL' TwinLabel ' - m_\nu^2 = 0 eV^2'];
katrinB  = ['KATRIN KSN1 ' CLflag '%CL' TwinLabel ' - free m_\nu^2 - 4 eV^2 pull term'];
katrinC  = ['KATRIN KSN1 ' CLflag '%CL' TwinLabel ' - m_\nu^2 = -1 eV^2'];
giunti   = 'arXiv:1912.12956 - 90%CL - [E_0-40;E_0+50] eV';

switch MainzTroitzk
    case 'ON'
        switch RealTwinFlag
            case 'Real'
                legend([p_m p_t pA1 pA2 pB pC],{mainz,troitsk,katrinA1,katrinA2,katrinB,katrinC},'Location','southwest','box','off');
            case 'Twin'
                legend([p_m p_t pA1 pA2 pB],{mainz,troitsk,katrinA1,katrinA2,katrinB},'Location','southwest','box','off')
        end
    case 'OFF'
        switch RealTwinFlag
            case 'Real'
                legend([pA1 pA2 pB pC],{katrinA1,katrinA2,katrinB,katrinC},'Location','southwest','box','off');
            case 'Twin'
                legend([pA1 pA2 pB],{katrinA1,katrinA2,katrinB},'Location','southwest','box','off')
        end
end
grid on
PRLFormat;

% Axis 
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');

axis([0.005 0.5 0.5 100])
axis square
%title('KATRIN Sterile Neutrino Analysis (KSN1) - 95% CL Sensitivity') % Exclusion Limit (Data)

fileString = sprintf('./plots/ksn1_exlimit_%sCL_m4sinSqT.pdf',CLflag);
export_fig(fig,fileString)