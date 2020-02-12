% configuration
TD = 'Flat60';
TimeSec =4*24*60*60;
MACE_Ba_T = 7*1e-04;
MACE_Bmax_T = 6*0.7;
WGTS_B_T    = 3.6*0.7;
WGTS_CD_MolPerCm2 = 5e17;

SysEffect = 'RF';
nTrials = 5000;
FSDNorm_RelErr = 0.01;
FSDShapeGS_RelErr = 0.04;
FSDShapeES_RelErr = 0.18;
MACE_Ba_T_RelErr = 0.002;
MACE_Bmax_T_RelErr = 0.002;
WGTS_B_T_RelErr = 0.002;
WGTS_CD_MolPerCm2_RelErr = 0.001;
ISXsection_RelErr = 0;
WGTS_TASR_RelErr = 0.001;

%% Init Model,  compute CM
A= ref_TBD_NominalKATRIN('TD',TD,'TimeSec',TimeSec,'recomputeRF','OFF',...
    'MACE_Ba_T',MACE_Ba_T,'WGTS_B_T',WGTS_B_T,'MACE_Bmax_T',MACE_Bmax_T,...
    'WGTS_CD_MolPerCm2',WGTS_CD_MolPerCm2);

[CMObj, MultiCM, MultiCMFrac, MultiCMShape, MultiCMNorm] = ...
    ref_CovarianceMatrix_NominalKATRIN('RecomputeFlag','OFF','ModelObj',A,...
    'nTrials',nTrials,'SysEffect',SysEffect,'PlotCM','OFF','SysBudget','99',...
    'FSDNorm_RelErr',FSDNorm_RelErr','FSDShapeGS_RelErr',FSDShapeGS_RelErr,...
    'FSDShapeES_RelErr',FSDShapeES_RelErr,'MACE_Ba_T_RelErr',MACE_Ba_T_RelErr,...
    'MACE_Bmax_T_RelErr',MACE_Bmax_T_RelErr,'WGTS_B_T_RelErr',WGTS_B_T_RelErr,...
    'WGTS_CD_MolPerCm2_RelErr',WGTS_CD_MolPerCm2_RelErr,'ISXsection_RelErr',ISXsection_RelErr,...
    'WGTS_TASR_RelErr',WGTS_TASR_RelErr);

switch SysEffect
    case 'RF'
CMObj.ComputeCM_RF;
    case 'FSD'
        CMObj.ComputeCM_FSD;
    case 'TASR'
        CMObj.ComputeCM_TASR;
    case 'TC'
        CMObj.ComputeCM_TC;
end
%% plot
CMObj.PlotCM('qUWindow',-5,'saveplot','OFF');
