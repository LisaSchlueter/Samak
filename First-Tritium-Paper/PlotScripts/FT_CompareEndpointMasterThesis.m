FT = MultiRunAnalysis('RunList','FTpaper','chi2','chi2Stat',...
    'ELossFlag','Abdurashitov','exclDataStart',13,'SysBudget',99,...%0
    'DataType','Real','FSDFlag','SAENZ',...
    'fixPar','1 5 6 7 8 9 10 11','NonPoissonScaleFactor',1,...
    'minuitOpt','min;migrad');

FT.ModelObj.DTFSD = 'DOSS';
FT.ModelObj.LoadFSD;
FT.ModelObj.ComputeTBDDS;
FT.ModelObj.ComputeTBDIS;

FT.ModelObj.WGTS_CD_MolPerCm2 = 4.56*1e17;
FT.ModelObj.ISXsection = 3.49e-22;
FT.ModelObj.AdjustRF;

FT.ComputeCM('BkgCM','OFF');

FT.Fit;