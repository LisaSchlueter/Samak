% ------------------------------------------------------------------------------------------------
% Goal: 
% Obtain sensitivity on neutrino mass squared
% Method: Monte Carlo Simulation
% Simulate KATRIN nSamples times (asimov + add stat fluctuations)
% Options: MTD, Background, Time
% Fit with Statistic or use a covariance matrix , 4 Parameter free(Nu,E0,N,B)
% Plot neutrino mass distribution
% Compute sensitivity to neutrino mass squared
% ------------------------------------------------------------------------------------------------
addpath(genpath('../../../Samak2.0'));

MACE_Ba_T_all = 7*1e-04;%(3:12)*1e-04;         % B-Field analyzing plane (T)
Range_all = 30;%[30, 45,60];                  % MTD energy range: 30,45,60eV below E0 (18575eV)    
nSamples = 1000;
TimeSec = 3*365*24*60*60;
SysEffect_all = {'FSD','RF','all','TC','TASR'};
chi2 = 'chi2CM';

if strcmp(chi2,'chi2Stat')
    SysEffect_all = '';
    SysFluct = 'OFF';
else
    SysFluct  = 'ON';
end
%for r=1:numel(Range_all)
range = Range_all;%Range_all(r);                 
%parfor i=1:numel(MACE_Ba_T_all)
MACE_Ba_T = MACE_Ba_T_all;%MACE_Ba_T_all(i);
TD = sprintf('Sensitivity_%.0feV_Ba%.0fG',range,MACE_Ba_T*1e4);
%% Do Monte Carlo
%stat
MC_SensitivityStudy_NominalKATRIN('saveResults','ON','StatFluct','ON','SysFluct','OFF',...
    'TimeSec',TimeSec,'chi2','chi2Stat','SysEffect','','nSamples',nSamples,...
    'TD',TD,'MACE_Ba_T',MACE_Ba_T,'mnuSq_i_Fit',0);

% systematics
parfor i=1:numel(SysEffect_all)
    SysEffect = SysEffect_all{i};
    MC_SensitivityStudy_NominalKATRIN('saveResults','ON','StatFluct','ON','SysFluct',SysFluct,...
        'TimeSec',TimeSec,'chi2',chi2,'SysEffect',SysEffect,'nSamples',nSamples,...
        'TD',TD,'MACE_Ba_T',MACE_Ba_T,'mnuSq_i_Fit',0);
end
%end
%end
