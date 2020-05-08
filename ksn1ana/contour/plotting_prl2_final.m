%% Plots for KSN1 sensitivity contours
%% In m4 sin2theta Plane

RealTwinFlag = 'Real';
MainzTroitzk = 'ON';

%% Data Importation

% Each datafile contains three variables : 
%   sith4_X (mixing angle)
%   mnu4Sq    (sterile mass)
%   chi_Z   (chiSq values, if everything is OK this should be only 4.61)

filepath   = [getenv('SamakPath'),'ksn1ana/contour/contourmatfiles/'];
CLflag     = '95'; % Confidence Level

switch RealTwinFlag
    case 'Real'
        TwinLabel   = '';
        % zero fixed nu mass - stat
        file_A1     = 'SamakContour_Real_95eV_chi2CMShape_E0BkgNorm.mat';
        % zero fixed nu mass - stat+sys
        file_A2     = 'SamakContour_Real_95eV_chi2CMShape_E0BkgNorm.mat';
        % free nu mass - stat+sys
        file_B      = 'SamakContour_Real_95eV_chi2CMShape_mNuE0BkgNorm_pull12.mat';
        % minus one fixed nu mass - stat+sys
        file_C      = 'SamakContour_Real_95eV_chi2CMShape_E0BkgNorm.mat';
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
db   = importdata([filepath,file_B]);    % Data B
dc   = importdata([filepath,file_C]);    % Data C


%% sin(th4)
da12  = da1.sin2T4_contour_95;
da22  = da2.sin2T4_contour_95;
db2   = db.sin2T4_contour_95;
dc2   = dc.sin2T4_contour_95;

% Constant data
d_giunti   = importdata([filepath,'coord_Giunti.mat']); % KATRIN Data from Giunti
d_mainz    = importdata([filepath,'coord_mainz.mat']);
d_troitsk  = importdata([filepath,'coord_troitsk.mat']);


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
%p_g     =       plot (d_giunti.sith4_X, d_giunti.mnu4Sq_contour_95,...
%                    'color',[0.9290 0.6940 0.1250],'LineWidth',1);
                
hold on

% === Our Data ===
% pA1      =       plot (([da12(cuta1),1]),(([da1.mnu4Sq_contour_95(cuta1),da1.mnu4Sq_contour_95(na1)])),...
%     ':','color',prlA,'LineWidth',4);

pA2      =       plot (([da22(cuta2),1]),(([da2.mnu4Sq_contour_95(cuta2),da2.mnu4Sq_contour_95(na2)])),...
    '-','color',prlA,'LineWidth',4);

hold on

plot([0.014 0.014],[62.5 62.5].^2,'s','LineWidth',3,'Color',prlA);

pB      =       plot (([db2(cutb),1]),(([db.mnu4Sq_contour_95(cutb),db.mnu4Sq_contour_95(nb)])),...
    '--','color',prlC,'LineWidth',4);

plot([0.014 0.014],[53 53].^2,'d','LineWidth',3,'Color',prlC);


switch RealTwinFlag
    case 'RealKNM1'
        hold on
        pC      =       plot (smooth([dc2(cutc),1]),smooth(([dc.mnu4Sq_contour_95(cutc),dc.mnu4Sq_contour_95(nc)])),...
            '-.','color',prlC,'LineWidth',4);
end

switch MainzTroitzk
    case 'ON'
        
        %   Mainz & Troitsk
        mainzSinTheta2 = 0.5 * (1-sqrt(1-d_mainz.sith4_X));
        p_m     =       plot (smooth(mainzSinTheta2), smooth((d_mainz.m4_Y)),...
            '-.','color',rgb('Salmon'),'LineWidth',2);
        % hold on
        troitskSinTheta2 = 0.5 * (1-sqrt(1-d_troitsk.sith4_X));
        p_t     =      plot (smooth(troitskSinTheta2), smooth((d_troitsk.m4_Y)),...
            '-.','color',rgb('DarkSlateGrey'),'LineWidth',2);
end

%% Plot Parameters

% Axis
xlabel('|U_{e4}|^2');
ylabel('m^2_{4}  (eV)');

% Labels
mainz    = 'Mainz 90%CL - m_\nu^2 = 0 eV^2';
troitsk  = 'Troitsk 90%CL - m_\nu^2 = 0 eV^2';
katrinA1  = ['KATRIN KSN1 ' CLflag '%CL' TwinLabel ' - stat only - m_\nu^2 = 0 eV^2'];
katrinA2  = ['KATRIN KSN1 ' CLflag '%CL' TwinLabel ' - m_\nu^2 = 0 eV^2'];
katrinB  = ['KATRIN KSN1 ' CLflag '%CL' TwinLabel ' - m_\nu^2 = nuisance parameter'];
katrinC  = ['KATRIN KSN1 ' CLflag '%CL' TwinLabel ' - m_\nu^2 = -1 eV^2'];
giunti   = 'arXiv:1912.12956 - 90%CL - [E_0-40;E_0+50] eV';

switch MainzTroitzk
    case 'ON'
        switch RealTwinFlag
            case 'Real'
                legend([p_m p_t  pA2 pB ],{mainz,troitsk,katrinA2,katrinB},'Location','southwest','box','off');
            case 'Twin'
                legend([p_m p_t  pA2 pB],{mainz,troitsk,katrinA2,katrinB},'Location','southwest','box','off')
        end
    case 'OFF'
        switch RealTwinFlag
            case 'Real'
                legend([pA1 pA2 pB ],{katrinA1,katrinA2,katrinB},'Location','southwest','box','off');
            case 'Twin'
                legend([pA1 pA2 pB],{katrinA1,katrinA2,katrinB},'Location','southwest','box','off')
        end
end
grid on
PRLFormat;

% Axis 
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');

axis([0.005 0.5 0.5 10000])
axis square
%title('KATRIN Sterile Neutrino Analysis (KSN1) - 95% CL Sensitivity') % Exclusion Limit (Data)

fileString = sprintf('./plots/ksn1_exlimit_%sCL_m4sinSqT.pdf',CLflag);
export_fig(fig,fileString)