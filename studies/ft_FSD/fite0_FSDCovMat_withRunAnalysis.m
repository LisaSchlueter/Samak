% script to do fits to data with FSD covmat (or other CMs) for Stacked Pixel

% General settings
myRunList = {'StackCD100all','StackCD100up','StackCD100down'};
ringCutFlag = 'ex2';
DataEffCorr = 'OFF';%'ROI+PileUp';
exclDataStart = 7;
% Systematics settings
nTrials = 1000;
RecomputeFlag = 'OFF';
StackCM = 'ON';
InitNormFit = 'ON';
%mySysEffects = struct('FSD','ON');
%Fit settings
fixPar = '1 5 6';
chi2Flag = 'chi2CM';
FSDNorm_RelErr = 0.01;
FSDShapeGS_RelErr = 0.04;
FSDShapeES_RelErr = 0.18;
WGTS_CD_MolPerCm2_RelErr = 0.08;
WGTS_TASR_RelErr = 0.001;

mySysEffects = struct(...
    'RF_EL','ON',...  % Response Function(RF) EnergyLoss
    'RF_BF','ON',...  % RF B-Fields
    'RF_RX','ON',...  % Column Density, inel cross ection
    'FSD','ON',...
    'TASR','ON',...
    'TCoff_RAD','ON',...
    'TCoff_OTHER','ON',...
    'DOPoff','OFF');
            
% Initialize RunAnalysis Model
A = MultiRunAnalysis('RunList',myRunList{1},'AnaFlag','StackPixel','DataEffCorr',DataEffCorr,...
    'fixPar',fixPar,'exclDataStart',exclDataStart);
if ~strcmp(chi2Flag,'chi2Stat')
A.ComputeCM('SysEffects',mySysEffects,'nTrials',nTrials,...
    'RecomputeFlag',RecomputeFlag,'InitNormFit',InitNormFit,...
    'FSDShapeGS_RelErr',FSDShapeGS_RelErr,'FSDShapeES_RelErr',FSDShapeES_RelErr,...
    'FSDNorm_RelErr',FSDNorm_RelErr,'WGTS_CD_MolPerCm2_RelErr',WGTS_CD_MolPerCm2_RelErr,...
    'WGTS_TASR_RelErr',WGTS_TASR_RelErr,'StackCM',StackCM);
end
A.chi2 = chi2Flag;
A.Fit;
close all;
A.PlotFit('ResidualsFlag','Norm');
