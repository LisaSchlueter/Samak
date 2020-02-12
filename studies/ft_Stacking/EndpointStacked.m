RunList=0;
close all;
%% 
if RunList == 0
%[ n, ~ , RunList ] = GetRunList( '../../tritium-data/hdf5/','2f-fpd*.h5',1,char('2f-fpd0040649.h5','40769.h5','40770.h5','40771.h5','40772.h5'));
[ n, ~ , RunList ] = GetRunList( '../../tritium-data/hdf5/','2f-fpd00*.h5',1,'string');
RunList=RunList(RunList>40531);             % First Good Run FT
%RunList=RunList(RunList<40693);            % 40538-40693 List of 100%CD 3h runs
Remove1=RunList==40773;RunList(Remove1)=[]; % Sterile Neutrino Run - High Count Rate
Remove2=RunList==40806;RunList(Remove2)=[]; % Sterile Neutrino Run - High Count Rate
Remove3=RunList==40995;RunList(Remove3)=[]; % Sterile Neutrino Run - High Count Rate
Remove4=RunList==40761;RunList(Remove4)=[]; % Stability Run - No Scan
Remove5=RunList==40762;RunList(Remove5)=[]; % Stability Run - No Scan
Remove6=RunList==40807;RunList(Remove6)=[]; % Anomaleous DT uncertainty
Remove7=RunList==40769;RunList(Remove7)=[]; % Sterile Neutrinos? 31 qU
Remove7=RunList==40770;RunList(Remove7)=[]; % Sterile Neutrinos? 31 qU
Remove8=RunList==40805;RunList(Remove8)=[]; % Sterile Neutrinos? 18 qU
RunList=RunList(RunList<40693);
%RunList=RunList(RunList<41003);
end
RunList = [40538:40543,40603,40604,40610:40613,40667:40693];

disp(RunList);
MRA = MultiRunAnalysis('RunList',RunList,'fixPar','1 5 6','chi2','chi2Stat',...
    'exclDataStart',9);

% Select Systematic Effects used in Fit
myEffects = struct(...
    'RF_EL','ON',...  % Response Function(RF) EnergyLoss
    'RF_BF','ON',...  % RF B-Fields
    'RF_RX','ON',...  % Column Density, inel cross ection
    'FSD','OFF',...
    'TASR','ON',...
    'TCoff_RAD','OFF',...
    'TCoff_OTHER','ON');
MRA.InitializeCM('SysEffects',myEffects,...
    'DataDriven','ON',...
    'WGTS_CD_MolPerCm2_RelErr',0.05);
E0    = zeros(7,1);
E0Err = zeros(7,1);
MRA.SimulateStackRuns;  %reset to init ModelObject
MRA.ComputeCM('InitNormFit','ON');
MRA.Fit;
MRA.PlotFit();

MRA.ComputeCM_Stack('Mode','Read');