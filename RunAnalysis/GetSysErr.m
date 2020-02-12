function [SysErr,CMArg] = GetSysErr(SysBudget)
% default systematic uncertainties, used in RunAnalysis and RunSensitivity
% do not change settings!!
% add new SysBudget number instead

if SysBudget==0 % default First Tritium 
    SysErr.WGTS_TASR_RelErr = 0.005;
    SysErr.FSDNorm_RelErr=  0.01;
    SysErr.FSDShapeGS_RelErr= 0.04;
    SysErr.FSDShapeES_RelErr= 0.18;
    SysErr.MACE_Ba_T_RelErr= 0.01;
    SysErr.MACE_Bmax_T_RelErr= 0.002;
    SysErr.WGTS_B_T_RelErr= 0.025;
    SysErr.WGTS_CD_MolPerCm2_RelErr= 0.03;
    SysErr.ISXsection_RelErr= 0.02;
    SysErr.DataDriven = 'ON';
    SysErr.FPDeff_RelErr = 0.001;
elseif SysBudget==99 % old First Tritium
    SysErr.WGTS_TASR_RelErr = 0.005;
    SysErr.FSDNorm_RelErr=  0.01;
    SysErr.FSDShapeGS_RelErr= 0.04;
    SysErr.FSDShapeES_RelErr= 0.18;
    SysErr.MACE_Ba_T_RelErr= 0.01;
    SysErr.MACE_Bmax_T_RelErr= 0.02;
    SysErr.WGTS_B_T_RelErr= 0.02;
    SysErr.WGTS_CD_MolPerCm2_RelErr= 0.022;
    SysErr.ISXsection_RelErr= 0.02;
    SysErr.DataDriven = 'OFF';
    SysErr.FPDeff_RelErr = 0.001;
elseif SysBudget==1
    SysErr.WGTS_TASR_RelErr = 5e-4;
    SysErr.FSDNorm_RelErr=  0.01;
    SysErr.FSDShapeGS_RelErr= 0.04;
    SysErr.FSDShapeES_RelErr= 0.18;
    SysErr.MACE_Ba_T_RelErr= 0.01;
    SysErr.MACE_Bmax_T_RelErr= 0.002;
    SysErr.WGTS_B_T_RelErr= 0.025;
    SysErr.WGTS_CD_MolPerCm2_RelErr= 0.0225;
    SysErr.ISXsection_RelErr= 0.006;
    SysErr.DataDriven = 'ON';
    SysErr.FPDeff_RelErr = 0.001;
elseif SysBudget ==2
    SysErr.WGTS_TASR_RelErr = 5e-4;
    SysErr.FSDNorm_RelErr=  0.01;
    SysErr.FSDShapeGS_RelErr= 0.04;
    SysErr.FSDShapeES_RelErr= 0.18;
    SysErr.MACE_Ba_T_RelErr= 0.01;
    SysErr.MACE_Bmax_T_RelErr= 0.002;
    SysErr.WGTS_B_T_RelErr= 0.025;
    SysErr.WGTS_CD_MolPerCm2_RelErr= 0.01;
    SysErr.ISXsection_RelErr= 0; %use rhod sigma together as uncertainty
    SysErr.DataDriven = 'ON';
    SysErr.FPDeff_RelErr = 0.001;
elseif SysBudget ==21
    SysErr.WGTS_TASR_RelErr = 5e-4;
    SysErr.FSDNorm_RelErr=  0.01;
    SysErr.FSDShapeGS_RelErr= 0.04;
    SysErr.FSDShapeES_RelErr= 0.18;
    SysErr.MACE_Ba_T_RelErr= 0.01;
    SysErr.MACE_Bmax_T_RelErr= 0.002;
    SysErr.WGTS_B_T_RelErr= 0.025;
    SysErr.WGTS_CD_MolPerCm2_RelErr= 0.008;
    SysErr.ISXsection_RelErr= 0; %use rhod sigma together as uncertainty
    SysErr.DataDriven = 'ON';
    SysErr.FPDeff_RelErr = 0.001;
elseif SysBudget ==22 % default Knm1
    SysErr.WGTS_TASR_RelErr = 5e-4;
    SysErr.FSDNorm_RelErr=  0.01;
    SysErr.FSDShapeGS_RelErr= 0.04;
    SysErr.FSDShapeES_RelErr= 0.18;
    SysErr.MACE_Ba_T_RelErr= 0.01;
    SysErr.MACE_Bmax_T_RelErr= 0.002;
    SysErr.WGTS_B_T_RelErr= 0.025;
    SysErr.WGTS_CD_MolPerCm2_RelErr= 0.0085;
    SysErr.ISXsection_RelErr= 0; %use rhod sigma together as uncertainty
    SysErr.DataDriven = 'ON';
    SysErr.FPDeff_RelErr = 1e-4;
elseif SysBudget == 31 % preliminary KNM2 systemaitcs (January 20)
    SysErr.WGTS_TASR_RelErr = 5e-4; % data driven
    SysErr.FSDNorm_RelErr=  0.01;
    SysErr.FSDShapeGS_RelErr= 0.04;
    SysErr.FSDShapeES_RelErr= 0.18;
    SysErr.MACE_Ba_T_RelErr= 0.004;
    SysErr.MACE_Bmax_T_RelErr= 0.002;
    SysErr.WGTS_B_T_RelErr= 0.025;
    SysErr.WGTS_CD_MolPerCm2_RelErr= 0.003;
    SysErr.ISXsection_RelErr= 0; %use rhod sigma together as uncertainty
    SysErr.DataDriven = 'ON';
    SysErr.FPDeff_RelErr = 1e-4;
end

CMArg = {'WGTS_CD_MolPerCm2_RelErr',SysErr.WGTS_CD_MolPerCm2_RelErr,...
    'MACE_Bmax_T_RelErr',SysErr.MACE_Bmax_T_RelErr,...
    'MACE_Ba_T_RelErr',SysErr.MACE_Ba_T_RelErr,...
    'WGTS_B_T_RelErr',SysErr.WGTS_B_T_RelErr,...
    'ISXsection_RelErr',SysErr.ISXsection_RelErr,...
    'WGTS_TASR_RelErr',SysErr.WGTS_TASR_RelErr,...
    'FSDNorm_RelErr',SysErr.FSDNorm_RelErr,...
    'FSDShapeGS_RelErr',SysErr.FSDShapeGS_RelErr,...
    'FSDShapeES_RelErr',SysErr.FSDShapeES_RelErr,...
    'FPDeff_RelErr',SysErr.FPDeff_RelErr,...
    'DataDriven',SysErr.DataDriven};
end
