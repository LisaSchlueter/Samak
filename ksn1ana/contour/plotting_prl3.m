%% KSN1 Letter Article 2020
%% Plots for KSN1 sensitivity contours
%% In Oscillation parameters Deltam41^2 sin^2(2theta) Plane
%% L.Schlueter / N. Le Guennic / T. Lasserre / 

% Range
Range      = '40';

% Confidence Level
CLflag     = '95'; % Confidence Level

% Path & KATRIN Results
filepath   = [getenv('SamakPath'),'ksn1ana/contour/contourmatfiles/'];
switch Range
    case '65'
        file_A     = 'SamakContour_Real_65eV_chi2CMShape_E0BkgNorm.mat';
    case '40'
        file_A    = 'SamakContour_Real_40eV_chi2CMShape_E0BkgNorm.mat';
end

da         = importdata([filepath,file_A]);   

% Conversion to Oscillation parameters Deltam41^2 sin^2(2theta)
switch CLflag
    case '95'
        ksn1SinSq2Theta  = 4*da.sin2T4_contour_95.*(1-da.sin2T4_contour_95);
        ksn1DeltaMSq     = da.mnu4Sq_contour_95;
    case '90'
        
        ksn1SinSq2Theta  = 4*da.sin2T4_contour_90.*(1-da.sin2T4_contour_90);
        ksn1DeltaMSq     = da.mnu4Sq_contour_90;
end

% Giunti's Results
d_giunti   = importdata([filepath,'coord_Giunti.mat']);        

% RAA
d_raa_90   = importdata([filepath,'coord_RAA_90.mat']);         
d_raa_95_a = importdata([filepath,'coord_RAA_95_A.mat']);
d_raa_95_b = importdata([filepath,'coord_RAA_95_B.mat']);

raa90x = d_raa_90.sith4_X;   raa90y = d_raa_90.m4_Y;   n90 = length(raa90x);
raa95x = d_raa_95_a.sith4_X; raa95y = d_raa_95_a.m4_Y; n95 = length(raa95x);

raa90x = [raa90x(1),raa90x,raa90x(n90)];
raa95x = [raa95x(1),raa95x,raa95x(n95)];
raa95y = [10000,raa95y,10000]; raa90y = [10000,raa90y,10000];

% Mainz / Troistk
switch CLflag
    case '95'
        d_mainz    = importdata([filepath,'coord_Mainz_95CL.mat']);
        d_troitsk  = importdata([filepath,'coord_Troitsk_95CL.mat']);
    case '90'
        d_mainz    = importdata([filepath,'coord_mainz.mat']);
        d_troitsk  = importdata([filepath,'coord_troitsk.mat']);
end

% Reactor Oscillation Experiments
d_stereo    = importdata([filepath,'coord_Stereo_95CL.mat']);
d_prospect  = importdata([filepath,'coord_Prospect_95CL.mat']);
d_neutrino4 = importdata([filepath,'coord_Neutrino4_123sigma.mat']);
d_danss     = importdata([filepath,'coord_DANSS_95CL.mat']);

%% Colors
prlG = [81 126 102]/255;
prlB = [50 148 216]/255;

%% Plotting
fig = figure('Renderer','painters');
set(fig, 'Units', 'normalized', 'Position', [0.001, 0.001,0.5, 0.999]);

%   Mainz & Troitsk
p_m     =       plot ((d_mainz.SinSquare2Theta_X), (d_mainz.DmSquare41_Y),...
                    '-.','color',rgb('Salmon'),'LineWidth',2);
hold on
p_t     =      plot ((d_troitsk.SinSquare2Theta_X), (d_troitsk.DmSquare41_Y),...
                    '--','color',rgb('DarkSlateGrey'),'LineWidth',2);
%   RAA
hold on
p_raa95 =       plot (raa95x,raa95y,...
                    'color',rgb('ForestGreen')*1.2,'LineWidth',2);
hold on
p_raa2  =       plot (d_raa_95_b.sith4_X, d_raa_95_b.m4_Y,...
                    'color',rgb('ForestGreen')*1.2,'LineWidth',2);
%   Giunti
%hold on
%p_g     =       plot (d_giunti.sith4_X, d_giunti.m4_Y,...
%                    'color',[0.9290 0.6940 0.1250],'LineWidth',1);
                
% Stereo 
hold on

p_s     =       plot (smooth(d_stereo.SinSquare2Theta_X),smooth(d_stereo.DmSquare41_Y),...
                    'color',rgb('GoldenRod'),'LineWidth',1);
                
% Prospect 
p_p     =       plot (smooth(d_prospect.SinSquare2Theta_X),smooth(d_prospect.DmSquare41_Y),...
                    'color',rgb('Violet'),'LineWidth',1);
                
% Neutrino4 
bn4     =       boundary((d_neutrino4.SinSquare2Theta_X_2sigma)',(d_neutrino4.DmSquare41_Y_2sigma)',0);
p_n     =       plot((d_neutrino4.SinSquare2Theta_X_2sigma(bn4)),(d_neutrino4.DmSquare41_Y_2sigma(bn4)),...
                    'color',rgb('Red'),'LineWidth',1);

% DANSS 
p_d     =       plot (smooth(d_danss.SinSquare2Theta_X),smooth(d_danss.DmSquare41_Y),...
                    'color',rgb('Grey'),'LineWidth',1);
                
% === KATRIN Results ===
hold on
pKNS1      =       plot (ksn1SinSq2Theta,ksn1DeltaMSq,'-','color',prlB,'LineWidth',5);

%% Plot Parameters

% Axis
xlabel('sin^2(2\theta_{ee})');
ylabel('\Deltam_{41}^2  (eV^2)');

% Labels
mainz     = ['Mainz ' CLflag '% CL'];
troitsk   = ['Troitsk ' CLflag '% CL'];
stereo    = 'Stéréo 95% CL';
prospect  = 'Prospect 95% CL';
neutrino4 = 'Neutrino-4 2\sigma';
danss     = 'DANSS 95% CL';
KatrinKSN1   = ['KATRIN KSN1 - stat+syst - ' CLflag '% CL'];% - [E_0-90;E_0+50] eV'];
giunti    = 'arXiv:1912.12956 - 90% CL - [E_0-40;E_0+50] eV';
raa90     = 'RAA+GA - PRD 83, 073006 (2011) - 90% CL';
raa95     = 'RAA+GA - PRD 83, 073006 (2011) - 95% CL';

hl=legend([p_m p_t p_s p_p p_n p_d pKNS1 p_raa95],...  % Label order omit p_g
        {mainz,troitsk,stereo,prospect,neutrino4,danss,...
        KatrinKSN1,...                             % KATRIN
        raa95},...                              % RAA
        'Location','northoutside',...           % Legend settings
        'box','off');
%        giunti,...                              % Giunti
hl.NumColumns=2;
grid on

PRLFormat;

% Axis 
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
axis([0.02 1 0.1 2000])
axis square

% Save Files
fileString = sprintf('./plots/ksn1_exlimit_%sCL_Dm24sinSq2T_%seVrange.pdf',CLflag,Range);
export_fig(fig,fileString)