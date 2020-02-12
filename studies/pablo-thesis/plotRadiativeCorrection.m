clear; close all;

% Add path to samak folder
addpath(genpath('../../../Samak2.0'));


% Corrections to the allowed diff. beta spectrum

% obj.TBDDS = obj.TBDDS .* obj.ComputeRadiativeCorr();
% obj.TBDDS = obj.TBDDS .* obj.ComputeRecoilWmVmACorr() ;
% obj.TBDDS = obj.TBDDS .* obj.ComputeFiniteExtChargeCorr() ;
% obj.TBDDS = obj.TBDDS .* obj.ComputeEEexchangeCorr(2) ;
% obj.TBDDS = obj.TBDDS .* obj.ComputeScreeningCorr() ;
% obj.TBDDS = obj.TBDDS .* obj.ComputeWintFiniteSizeCorr() ;
% obj.TBDDS = obj.TBDDS .* obj.ComputeRecoilCoulombCorr() ;
 


A = ref_corrections();
%A.ComputeTBDDS();

B = ref_corrections('RadiativeFlag','ON','qUmin',18575-1600,'qUmax',18574,'Mode','Sim',...
    'nTeBinningFactor',1,'KTFFlag','OFF');
B.Q_i = 18575;
B.SetFitBias(0)
B.SetKinVariables();
B.RadType = 1;
%B.ComputeTBDDS();

% relativeDiff = (B.TBDDS)./(A.TBDDS);
relativeDiff = B.ComputeRadiativeCorr();

plot(B.Te(B.Te<18575)-B.Q_i,relativeDiff(B.Te<18575),'k','LineWidth',3)

title('Radiative Correction (G)');
xlabel('E - E_0 (eV)');
ylabel('rel. change in spectrum S_G/S');

max_change_FT = max(relativeDiff(B.Te<18575))-min(relativeDiff(B.Te<18575));
max_change_60 = max(relativeDiff(B.Te>18575-60))-min(relativeDiff(B.Te>18575-60));

legend_string = sprintf('max. change in FT range =  %.g %% \nmax. change -60 eV below E_0 = %.g %%',...
    max_change_FT*100,max_change_60*100);
legend(legend_string,'location','southwest')

PrettyFigureFormat;

export_fig('plots/RadiativeCorrections.pdf','-pdf')

