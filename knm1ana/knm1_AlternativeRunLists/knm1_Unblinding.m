
try
    knm1_RandomHalf;
catch
end

try
    ColumnDensityScan('chi2Stat')
catch
end
try
    ColumnDensityScan('chi2CMShape')
catch
end

try
    M = MultiRunAnalysis('RunList','KNM1','DataType','Twin','exclDataStart',14,'fixPar','5 6 7 8 9 10 11',...
        'SysBudget',22,'chi2','chi2CMShape');
    S = RunSensitivity('RunAnaObj',M);
    [UpperLimit,mNuSq_Quantil90,mNuSq_t]  = ...
        S.ComputeUpperLimit('nSamples',1000,'mNuSq_t',(0.95)^2);%[0,0.5,1.2]1.2:-0.2:0.2);%[0,0.5,1,1.5,2,3,0.2:0.2:1.2]);
catch
    
end

try
    M = MultiRunAnalysis('RunList','KNM1','DataType','Twin','exclDataStart',14,'fixPar','5 6 7 8 9 10 11',...
        'SysBudget',22,'chi2','chi2CMShape','TwinBias_mnuSq',(0.94)^2);
    S = RunSensitivity('RunAnaObj',M);
    S.PlotCompareSensitivityTwinData;
catch
end
