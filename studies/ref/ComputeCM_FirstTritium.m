%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to be run on server
% compute covariance matrix for all runs individually with exact slow
% control data and exact column density
addpath(genpath('../../../Samak2.0'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model Options
myRunNr = GetFTRunList('StackCD100all');
myringcutFlag= 'ex2';        % exclude 2 outer rings
AnaFlag = 'StackPixel';
myDataEffCorr ='ROI+PileUp'; % PixelWise Efficiency Correction of Data
% Covariance Matrix Options
mySysEffects =  struct(...
                'RF_EL','ON',...  % Response Function(RF) EnergyLoss
                'RF_BF','ON',...  % RF B-Fields
                'RF_RX','ON',...  % Column Density, inel cross ection
                'FSD','OFF',...
                'TASR','OFF',...
                'TCoff_RAD','OFF',...
                'TCoff_OTHER','OFF',...
                'DOPoff','OFF');
nTrials = 1000;            
RecomputeFlag = 'OFF';
DataDriven = 'OFF'; 
WGTS_CD_MolPerCm2_RelErr= 0.08;

for i=1:numel(myRunNr)
AnaObj = RunAnalysis(...
    'RunNr',myRunNr(i),'ringCutFlag',myringcutFlag,...
    'DataEffCorr',myDataEffCorr,'AnaFlag',AnaFlag);
AnaObj.ComputeCM(...
    'SysEffects',mySysEffects,'nTrials',nTrials,'RecomputeFlag',RecomputeFlag,...
    'DataDriven',DataDriven,'WGTS_CD_MolPerCm2_RelErr',WGTS_CD_MolPerCm2_RelErr);
%AnaObj.FitCM_Obj.PlotCM('saveplot','OFF');
end