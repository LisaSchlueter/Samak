%
% WGTS Potential Radial Profile
% Suggested by Manuel Klein
% 29/10/2019
%

% Inputs
WGTS_Radius_cm = 4.5;    % WGTS Tube Radius in cm
U_RW           = 1; % RW Potential  
SPR1_cm        = (1/4-1/8);
SPR2_cm        = (2/4-1/8).*WGTS_Radius_cm;
SPR3_cm        = (3/4-1/8).*WGTS_Radius_cm;
SPR4_cm        = (4/4-1/8).*WGTS_Radius_cm;

%  WGTS Potential Radial Profile
% U_WGTS   = @(r) U_RW .* (1 - (r./WGTS_Radius_cm).^2);
U_WGTS   = @(r,offset) U_RW .* (1 - (1-offset)*(r./WGTS_Radius_cm).^2);

%  WGTS Potential Radial Profile
U_Offset = @(coeff) coeff .* (U_RW - U_WGTS([SPR1_cm SPR2_cm SPR3_cm SPR4_cm],0));

%% Display
radius = 0:0.1:WGTS_Radius_cm;
f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
plot(radius,U_WGTS(radius,0)*100,'LineWidth',5);
%plot(1:4,U_Offset(1),'LineWidth',5);
xlabel('Radius of WGTS Beam Tube');
ylabel(sprintf('Rel. plasma potential (%%)'));
%ylabel('WGTS Plasma Potential (V)');
xlim([0,4.5]);
PrettyFigureFormat('FontSize',26);
export_fig(f1,[getenv('SamakPath'),'knm2ana/knm2_RingFit/plots/WGTS_PotentialRadialProfile.pdf']);


