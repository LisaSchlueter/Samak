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
 


A = ref_corrections('qUmin',18575-1600,'qUmax',18570,'Mode','Sim',...
    'nTeBinningFactor',1,'KTFFlag','OFF');
A.Q_i = 18575;
A.SetFitBias(0)
A.SetKinVariables()
A.ComputeTBDDS();

AllCorr = {'ScreeningFlag','OFF','FiniteExtChargeFlag','ON',...
'WintFiniteSizeFlag','ON','EEexchangeFlag','ON','RecoilCoulombFlag','ON',...
'RecoilWmVmAFlag','ON'};
            

B = ref_corrections('RadiativeFlag','OFF','qUmin',18575-1600,'qUmax',18570,'Mode','Sim',...
    'nTeBinningFactor',1,'KTFFlag','OFF',AllCorr{:});
B.Q_i = 18575;
B.SetFitBias(0)
B.SetKinVariables();
B.RadType = 1;
B.ComputeTBDDS();

relativeDiff = (A.TBDDS)./(B.TBDDS);
%relativeDiff = B.ComputeRadiativeCorr();

plot(B.Te(B.Te<18575)-18575,relativeDiff(B.Te<18575),'b','LineWidth',3)

title('All Corrections except radiative');
xlabel('E - E_0 (eV)');
ylabel('relative change in spectrum S_{All}/S')

PrettyFigureFormat;
max(relativeDiff(B.Te<18575))-min(relativeDiff(B.Te<18575))
max(relativeDiff(B.Te>18575-60))-min(relativeDiff(B.Te>18575-60))

export_fig('AllCorr.pdf','-pdf');
