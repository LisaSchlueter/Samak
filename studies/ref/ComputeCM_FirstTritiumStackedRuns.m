%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to be run on server
% compute covariance matrix for all runs individually with exact slow
% control data and exact column density
addpath(genpath('../../../Samak2.0'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model Options
%myList = [40539:40543,40603,40604,40610:40613,40667:40693,40976,40977,40979,40980,40982,40983,40985,40986,40988,40989,40991,40992,40994,40997,40998,41001,41002];%:41033];%,41045:41048];
%myList = {'StackCD100all','StackCD100up','StackCD100down'};
myList = 'StackCD100all_3hours';
myringcutFlag= 'ex2';        % exclude 2 outer rings
AnaFlag = 'StackPixel';
myDataEffCorr ='ROI+PileUp'; % PixelWise Efficiency Correction of Data
% Covariance Matrix Options
mySysEffects =  struct(...
                'RF_EL','ON',...  % Response Function(RF) EnergyLoss
                'RF_BF','ON',...  % RF B-Fields
                'RF_RX','ON',...  % Column Density, inel cross ection
                'FSD','ON',...
                'TASR','ON',...
                'TCoff_RAD','ON',...
                'TCoff_OTHER','ON',...
                'DOPoff','OFF');
nTrials = 1000;            
RecomputeFlag = 'OFF';
DataDriven = 'OFF'; 
WGTS_CD_MolPerCm2_RelErr= 0.08;
FSDShapeGS_RelErr = 0.04;
FSDShapeES_RelErr = 0.18;
WGTS_TASR_RelErr = 0.001;

%parfor i=1:numel(myList)
AnaObj = MultiRunAnalysis('RunList',myList,'ringCutFlag',myringcutFlag,'DataEffCorr',myDataEffCorr,...
    'AnaFlag',AnaFlag);
AnaObj.ComputeCM('SysEffects',mySysEffects,'nTrials',nTrials,'RecomputeFlag',RecomputeFlag,...
   'DataDriven',DataDriven,'WGTS_CD_MolPerCm2_RelErr',WGTS_CD_MolPerCm2_RelErr,...
    'FSDShapeGS_RelErr',FSDShapeGS_RelErr,'FSDShapeES_RelErr',FSDShapeES_RelErr,...
    'WGTS_TASR_RelErr',WGTS_TASR_RelErr);
%AnaObj.FitCM_Obj.PlotStack;
%end
