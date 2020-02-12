% Computation of response function covariance matrix without energy loss
% function uncertainty
clear;
addpath(genpath('../../../Samak2.0'));

WGTS_CD_MolPerCm2_RelErr = 0.05;
ISXsection_RelErr = 0.02;
WGTS_B_T_RelErr = 0.02;
MACE_Bmax_T_RelErr = 0.02;
MACE_Ba_T_RelErr = 0.02;

myEffects = struct(...
    'RF_EL','OFF',...  % Response Function(RF) EnergyLoss
    'RF_BF','ON',...  % RF B-Fields
    'RF_RX','ON',...  % Column Density, inel cross ection
    'FSD','OFF',...
    'TASR','OFF',...
    'TCoff_RAD','OFF',...
    'TCoff_OTHER','OFF',...
    'DOPoff','OFF');
R = MultiRunAnalysis('RunList','StackCD100_3hours','fixPar','1 5 6');
R.ComputeCM('SysEffects',myEffects,'MACE_Bmax_T_RelErr',MACE_Bmax_T_RelErr,...
             'MACE_Ba_T_RelErr',MACE_Ba_T_RelErr,'MACE_Ba_T_RelErr',MACE_Ba_T_RelErr,...
             'ISXsection_RelErr',ISXsection_RelErr,'WGTS_CD_MolPerCm2_RelErr',WGTS_CD_MolPerCm2_RelErr);