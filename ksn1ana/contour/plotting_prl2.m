%% KSN1 Letter Article 2020
%% Plots for KSN1 sensitivity contours
%% In m4 sin2theta Plane
%% L.Schlueter / N. Le Guennic / T. Lasserre / 

% Range
Range        = '65';

% Plot Real = Data or Twin
RealTwinFlag = 'Real';

% Plot Mainz/Troistk
MainzTroitsk = 'ON';

% Confidence Level
CLflag       = '95'; 

%% Data Importation
filepath   = [getenv('SamakPath'),'ksn1ana/contour/contourmatfiles/'];
switch RealTwinFlag
    case 'Real'
        TwinLabel   = '';
        
        switch Range
            case '65'
        % KATRIN KSN1 zero fixed nu mass - stat - 65eV
        file_A1     = 'SamakContour_Real_65eV_chi2CMShape_E0BkgNorm.mat';
        
        % KATRIN KSN1 zero fixed nu mass - stat+sys - 65eV
        file_A2     = 'SamakContour_Real_65eV_chi2CMShape_E0BkgNorm.mat';
        BFsin2theta = 0.0277;
        BFm24       = 68.6182; 
        
        % KATRIN KSN1 free nu mass - stat+sys  - 65eV
        % file_B      = 'SamakContour_Real_65eV_chi2CMShape_mNuE0BkgNorm.mat';
        % BFsin2theta_mfree = 0.0182;
        % BFm24_mfree       = 206.5607; % eV^2
        
        % KATRIN KSN1 free nu mass - stat+sys + pull - 65eV
        file_B      = 'SamakContour_Real_65eV_chi2CMShape_mNuE0BkgNorm_pull12.mat';
        BFsin2theta_mfree = 0.022;
        BFm24_mfree       = 115.7005; % eV^2
        
            case '40'
                
        % KATRIN KSN1 zero fixed nu mass - stat - 40eV
        file_A1     = 'SamakContour_Real_40eV_chi2CMShape_E0BkgNorm.mat';
        
        % KATRIN KSN1 zero fixed nu mass - stat+sys - 40eV
        file_A2     = 'SamakContour_Real_40eV_chi2CMShape_E0BkgNorm.mat';
        BFsin2theta = 0;
        BFm24       = 0; 
        
        % KATRIN KSN1 free nu mass - stat+sys + pull - 40eV
        file_B      = 'SamakContour_Real_40eV_chi2CMShape_E0BkgNorm.mat';
        BFsin2theta_mfree = 0;
        BFm24_mfree       = 0; % eV^2
        
        end
        
    case 'Twin'
        TwinLabel   = 'Simulation';
        
        % KATRIN KSN1 zero fixed nu mass - stat
        file_A1     = 'SamakContour_Twin_65eV_chi2CMShape_E0BkgNorm.mat';
        
        % KATRIN KSN1 zero fixed nu mass - stat+sys
        file_A2     = 'SamakContour_Twin_65eV_chi2CMShape_E0BkgNorm.mat';
        BFsin2theta = 0.;
        BFm24       = 0.; 
        
        % KATRIN KSN1 free nu mass - stat+sys 
        % file_B      = 'SamakContour_Twin_65eV_chi2CMShape_mNuE0BkgNorm.mat';
        % BFsin2theta_mfree = 0.;
        % BFm24_mfree       = 0.; % eV^2
        
        % KATRIN KSN1 free nu mass - stat+sys + pull
        file_B      = 'SamakContour_Twin_65eV_chi2CMShape_mNuE0BkgNorm_pull12.mat';
        BFsin2theta_mfree = 0.;
        BFm24_mfree       = 0.; % eV^2
end

da1  = importdata([filepath,file_A1]);   % Data A1
da2  = importdata([filepath,file_A2]);   % Data A2
db   = importdata([filepath,file_B]);    % Data B

% sin^2(theta14) = |Ue4|^2
switch CLflag
    case '95'
da12  = da1.sin2T4_contour_95;
da22  = da2.sin2T4_contour_95;
db2   = db.sin2T4_contour_95;
    case '90'
da12  = da1.sin2T4_contour_90;
da22  = da2.sin2T4_contour_90;
db2   = db.sin2T4_contour_90;
end


% Giunti Results
d_giunti   = importdata([filepath,'coord_Giunti.mat']); 

% Mainz / Troistk
d_mainz    = importdata([filepath,'coord_Mainz_95CL.mat']);
d_troitsk  = importdata([filepath,'coord_Troitsk_95CL.mat']);

% Plot tunings
na1   = length(da12); na2   = length(da22); nb   = length(db2);  %nc   = length(dc2);  
cuta1 = (1:na1); cuta2 = (1:na2);  cutb = (1:nb);  %cutc = (1:nc);  

%% Colors
prlA = [50 148 216]/255;
prlB = rgb('IndianRed');
prlC = rgb('DarkGreen');

%% Plotting
fig = figure('Renderer','painters'); 
set(fig, 'Units', 'normalized', 'Position', [0.001, 0.001,0.5, 0.999]);

%hold on
%   Giunti
%p_g     =       plot (d_giunti.sith4_X, d_giunti.mnu4Sq_contour_95,...
%                    'color',[0.9290 0.6940 0.1250],'LineWidth',1);
                
hold on

% === KATRIN Results ===

% pA1      =       plot (([da12(cuta1),1]),(([da1.mnu4Sq_contour_95(cuta1),da1.mnu4Sq_contour_95(na1)])),...
%     ':','color',prlA,'LineWidth',4);

pA2      =       plot (([da22(cuta2),1]),(([da2.mnu4Sq_contour_95(cuta2),da2.mnu4Sq_contour_95(na2)])),...
    '-','color',prlA,'LineWidth',4);

hold on

%plot(BFsin2theta,BFm24,'s','LineWidth',3,'Color',prlA);

pB      =       plot (([db2(cutb),1]),(([db.mnu4Sq_contour_95(cutb),db.mnu4Sq_contour_95(nb)])),...
    ':','color',prlC,'LineWidth',4);

%plot(BFsin2theta_mfree,BFm24_mfree,'+','LineWidth',3,'Color',prlC);

% === Mainz / Troitsk Results ===

switch MainzTroitsk
    case 'ON'
        %   Mainz & Troitsk
        mainzSinTheta2 = 0.5 * (1-sqrt(1-d_mainz.SinSquare2Theta_X));
        p_m     =       plot ((mainzSinTheta2), ((d_mainz.DmSquare41_Y)),...
            '-.','color',rgb('Salmon'),'LineWidth',2);
        % hold on
        troitskSinTheta2 = 0.5 * (1-sqrt(1-d_troitsk.SinSquare2Theta_X));
        p_t     =      plot ((troitskSinTheta2), ((d_troitsk.DmSquare41_Y)),...
            '-.','color',rgb('DarkSlateGrey'),'LineWidth',2);
end

%% Plot Parameters

% Axis
xlabel('|U_{e4}|^2');
ylabel('m^2_{4}  (eV)');

% Labels
mainz     = 'Mainz 95%CL - m_\nu^2 = 0 eV^2';
troitsk   = 'Troitsk 95%CL - m_\nu^2 = 0 eV^2';
katrinA1  = [TwinLabel ' KATRIN KSN1 ' CLflag '%CL - stat only - m_\nu^2 = 0 eV^2'];
katrinA2  = [TwinLabel ' KATRIN KSN1 ' CLflag '%CL - m_\nu^2 = 0 eV^2'];
katrinB   = [TwinLabel ' KATRIN KSN1 ' CLflag '%CL - m_\nu^2 free - \sigma(m_\nu^2)=1.94 eV^2'];
giunti    = 'arXiv:1912.12956 - 90%CL - [E_0-40;E_0+50] eV';

switch MainzTroitsk
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

axis([0.005 0.5 1 5000])
axis square

fileString = sprintf('./plots/%sksn1_exlimit_%sCL_m4sinSqT_%seVrange.pdf',TwinLabel,CLflag,Range);
export_fig(fig,fileString)