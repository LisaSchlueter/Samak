function SysBudgetInfo = GetSysBudget(varargin)
p = inputParser;
p.addParameter('SysBudget','05',@(x)ischar(x));
p.addParameter('ELossFlag','Aseev',@(x)ischar(x));
p.parse(varargin{:});
SysBudget = p.Results.SysBudget;
ELossFlag = p.Results.ELossFlag;

switch SysBudget
    case '00'    
        SysBudgetInfo.MACE_Ba_T_RelErr = 0.002;
        SysBudgetInfo.MACE_Ba_T_Err_str = sprintf('0.2%%');
        SysBudgetInfo.MACE_Bmax_T_RelErr = 0.002;
        SysBudgetInfo.WGTS_B_T_RelErr = 0.002;       
        SysBudgetInfo.WGTS_CD_MolPerCm2_RelErr = 0.002;
        SysBudgetInfo.ISXsection_RelErr = 0.02;    
        SysBudgetInfo.WGTS_TASR_RelErr = 0.001;
        SysBudgetInfo.FSDNorm_RelErr = 0.01;
        SysBudgetInfo.FSDShapeGS_RelErr = 0.04;
        SysBudgetInfo.FSDShapeES_RelErr = 0.18;  
    case '01'
        SysBudgetInfo.WGTS_B_T_RelErr    = 0.002;
        SysBudgetInfo.MACE_Ba_T_Err   = 2*1e-06; %absolute value!
        SysBudgetInfo.MACE_Ba_T_Err_str = sprintf('2 \\muT');
        SysBudgetInfo.MACE_Bmax_T_RelErr = 0.002;
        SysBudgetInfo.WGTS_CD_MolPerCm2_RelErr = 0.001;
        SysBudgetInfo.ISXsection_RelErr = 0.000;
        SysBudgetInfo.WGTS_TASR_RelErr = 0.000;
        SysBudgetInfo.FSDNorm_RelErr = 0.010;
        SysBudgetInfo.FSDShapeGS_RelErr = 0.040;
        SysBudgetInfo.FSDShapeES_RelErr = 0.180;
    case '02'
        SysBudgetInfo.WGTS_B_T_RelErr    = 0.002;
        SysBudgetInfo.MACE_Ba_T_Err   = 2*1e-06; %absolute value!
        SysBudgetInfo.MACE_Ba_T_Err_str = sprintf('2 \\muT');
        SysBudgetInfo.MACE_Bmax_T_RelErr = 0.002;
        SysBudgetInfo.WGTS_CD_MolPerCm2_RelErr = 0.002;
        SysBudgetInfo.ISXsection_RelErr = 0.0046;
        SysBudgetInfo.WGTS_TASR_RelErr = 0.000;
        SysBudgetInfo.FSDNorm_RelErr = 0.010;
        SysBudgetInfo.FSDShapeGS_RelErr = 0.040;
        SysBudgetInfo.FSDShapeES_RelErr = 0.180;
    case '03'
        SysBudgetInfo.MACE_Ba_T_Err_str = sprintf('Ba\\cdot1.9\\cdot10^{-3}+2.6\\cdot10^{-6} T');
        SysBudgetInfo.MACE_Bmax_T_RelErr = 0.002;
        SysBudgetInfo.WGTS_B_T_RelErr = 0.002;
        SysBudgetInfo.WGTS_CD_MolPerCm2_RelErr = 0.001;
        SysBudgetInfo.ISXsection_RelErr = 0;
        SysBudgetInfo.WGTS_TASR_RelErr = 0;
        SysBudgetInfo.FSDNorm_RelErr = 0.01;
        SysBudgetInfo.FSDShapeGS_RelErr = 0.04;
        SysBudgetInfo.FSDShapeES_RelErr = 0.18;   
    case '04'
        SysBudgetInfo.MACE_Ba_T_Err_str = sprintf('Ba\\cdot1.9\\cdot10^{-3}+2.6\\cdot10^{-6} T');
        SysBudgetInfo.MACE_Bmax_T_RelErr = 0.002;
        SysBudgetInfo.WGTS_B_T_RelErr = 0.002;
        SysBudgetInfo.WGTS_CD_MolPerCm2_RelErr = 0.002;
        SysBudgetInfo.ISXsection_RelErr = 0.0046;
        SysBudgetInfo.WGTS_TASR_RelErr = 0;
        SysBudgetInfo.FSDNorm_RelErr = 0.01;
        SysBudgetInfo.FSDShapeGS_RelErr = 0.04;
        SysBudgetInfo.FSDShapeES_RelErr = 0.18;
    case '05'
        SysBudgetInfo.MACE_Ba_T_Err_str = sprintf('0.2%%');
        SysBudgetInfo.MACE_Bmax_T_RelErr = 0.002;
        SysBudgetInfo.WGTS_B_T_RelErr = 0.002;
        SysBudgetInfo.WGTS_CD_MolPerCm2_RelErr = 0.0014;
        SysBudgetInfo.ISXsection_RelErr = 0.0014;
        SysBudgetInfo.WGTS_TASR_RelErr = 0;
        SysBudgetInfo.FSDNorm_RelErr = 0.01;
        SysBudgetInfo.FSDShapeGS_RelErr = 0.04;
        SysBudgetInfo.FSDShapeES_RelErr = 0.18;
    case '06'% Thierry - KATRIN 2019 1 month
        SysBudgetInfo.MACE_Ba_T_Err_str = sprintf('Ba\\cdot1.9\\cdot10^{-3}+2.6\\cdot10^{-6} T');
        SysBudgetInfo.MACE_Bmax_T_RelErr = 0.01;
        SysBudgetInfo.WGTS_B_T_RelErr = 0.01;
        SysBudgetInfo.WGTS_CD_MolPerCm2_RelErr = 0.02;
        SysBudgetInfo.ISXsection_RelErr = 0.0046;
        SysBudgetInfo.WGTS_TASR_RelErr = 0.005;
        SysBudgetInfo.FSDNorm_RelErr = 0.01;
        SysBudgetInfo.FSDShapeGS_RelErr = 0.04;
        SysBudgetInfo.FSDShapeES_RelErr = 0.18;
    case '07'% Thierry - KATRIN 2021 1 year
        SysBudgetInfo.MACE_Ba_T_Err_str = sprintf('Ba\\cdot1.9\\cdot10^{-3}+2.6\\cdot10^{-6} T');
        SysBudgetInfo.MACE_Bmax_T_RelErr = 0.005;
        SysBudgetInfo.WGTS_B_T_RelErr = 0.005;
        SysBudgetInfo.WGTS_CD_MolPerCm2_RelErr = 0.005;
        SysBudgetInfo.ISXsection_RelErr = 0.0046;
        SysBudgetInfo.WGTS_TASR_RelErr = 0.005;
        SysBudgetInfo.FSDNorm_RelErr = 0.01;
        SysBudgetInfo.FSDShapeGS_RelErr = 0.04;
        SysBudgetInfo.FSDShapeES_RelErr = 0.18;
    case '08'% Thierry - KATRIN 2023 - No E-loss 3 years
        SysBudgetInfo.MACE_Ba_T_Err_str = sprintf('Ba\\cdot1.9\\cdot10^{-3}+2.6\\cdot10^{-6} T');
        SysBudgetInfo.MACE_Bmax_T_RelErr = 0.002;
        SysBudgetInfo.WGTS_B_T_RelErr = 0.002;
        SysBudgetInfo.WGTS_CD_MolPerCm2_RelErr = 0.002;
        SysBudgetInfo.ISXsection_RelErr = 0.0046;
        SysBudgetInfo.WGTS_TASR_RelErr = 0.005;
        SysBudgetInfo.FSDNorm_RelErr = 0.01;
        SysBudgetInfo.FSDShapeGS_RelErr = 0.04;
        SysBudgetInfo.FSDShapeES_RelErr = 0.18;
    case '09'% Thierry - KATRIN 2023 - No E-loss
        SysBudgetInfo.MACE_Ba_T_Err_str = sprintf('Ba\\cdot1.9\\cdot10^{-3}+2.6\\cdot10^{-6} T');
        SysBudgetInfo.MACE_Bmax_T_RelErr = 0.002;
        SysBudgetInfo.WGTS_B_T_RelErr = 0.002;
        SysBudgetInfo.WGTS_CD_MolPerCm2_RelErr = 0.001;
        SysBudgetInfo.ISXsection_RelErr = 0;
        SysBudgetInfo.WGTS_TASR_RelErr = 0;
        SysBudgetInfo.FSDNorm_RelErr = 0.01;
        SysBudgetInfo.FSDShapeGS_RelErr = 0.04;
        SysBudgetInfo.FSDShapeES_RelErr = 0.18;
              
end

if ~strcmp(SysBudget,'08') && ~strcmp(SysBudget,'09') % if sysbudget NOT 09 or 08
    if strcmp(ELossFlag,'Aseev')
    SysBudgetInfo.ELoss = sprintf('\nEnergy Loss Uncertainties from Eur. Phys. J. D 10, 39-52 (Aseev et.al.)'); 
    elseif strcmp(ELossFlag,'Abdurashitov')
     SysBudgetInfo.ELoss  =  sprintf('\nEnergy Loss Uncertainties from Phys.Part.Nucl.Lett. 14 (2017) no.6, 892-899 (Abdurashitov et al)');    
    end
elseif strcmp(SysBudget,'08') || strcmp(SysBudget,'09') % if sysbudget 09 or 08
    SysBudgetInfo.ELoss = sprintf('\nEnergy Loss Uncertainties negligible');
end
end