%% calculate grid for ksn1 stat. + 1 syst
range = 95;
nGridSteps = 25;
chi2Str = 'chi2CMShape';
DataType = 'Twin';
freePar = 'E0 Bkg Norm';
RunList = 'KNM1';
SmartGrid = 'OFF';

mySysEffects = {...
    'TASR',...             % tritium activity fluctuations from subrun to subrun
    'RF_EL',...            % energy-loss
    'RF_BF',...            % B-fields
    'RF_RX',...            % column Density, inel. cross ection
    'FSD',...              % final state distribution       
    'TCoff_OTHER',...      % theoretical corrections without radiative
    'Stack',...            % stacking
    'FPDeff'};             % fpd efficiecny

%% load grid (or calculate if doesn't exist)
for i=1:numel(mySysEffects)
    [mnu4Sq,sin2T4,chi2,chi2_ref,savefile] = KSN1GridSearch('range',range,...
        'nGridSteps',nGridSteps,...
        'chi2',chi2Str,...
        'DataType',DataType,...
        'freePar',freePar,...
        'RunList',RunList,...
        'SmartGrid',SmartGrid,...
        'RecomputeFlag','OFF',...
        'SysEffect',mySysEffects{i});
end