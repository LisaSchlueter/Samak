% neutrino mass shift due to magnetic fields
RunList = 'KNM1';
mNu = zeros(3,1); mNuErr = zeros(3,1);
if ~exist('A','var')
    A = MultiRunAnalysis('RunList',RunList,'chi2','chi2Stat','DataType','Twin','SysBudget',1,'exclDataStart',2,'fixPar','5 6 7 8 9 10 11');
    A.Fit;
    mNu(1)    = A.FitResult.par(1);
    mNuErr(1) = A.FitResult.err(1);
    MACE_Bmax_T_i = A.ModelObj.MACE_Bmax_T;
    MACE_Ba_T_i   = A.ModelObj.MACE_Ba_T;
    WGTS_B_T_i    = A.ModelObj.WGTS_B_T;
end

SysErr = GetSysErr(A.SysBudget);
A.SimulateStackRuns('MACE_Ba_T',MACE_Ba_T_i.*(1+SysErr.MACE_Ba_T_RelErr),...
                    'MACE_Bmax_T',MACE_Bmax_T_i.*(1+SysErr.MACE_Bmax_T_RelErr),...
                    'WGTS_B_T',WGTS_B_T_i.*(1+SysErr.WGTS_B_T_RelErr));
A.Fit;
mNu(2)   = A.FitResult.par(1);
mNuErr(2) = A.FitResult.err(1);

A.SimulateStackRuns('MACE_Ba_T',MACE_Ba_T_i.*(1-SysErr.MACE_Ba_T_RelErr),...
                    'MACE_Bmax_T',MACE_Bmax_T_i.*(1-SysErr.MACE_Bmax_T_RelErr),...
                    'WGTS_B_T',WGTS_B_T_i.*(1-SysErr.WGTS_B_T_RelErr));
A.Fit;
mNu(3)   = A.FitResult.par(1);
mNuErr(3) = A.FitResult.err(1);
%%
