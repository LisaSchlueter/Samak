MRA = MultiRunAnalysis('RunList','StackCD100all');
MRA.exclDataStart = 1;
MRA.ModelObj.ComputeTBDDS; MRA.ModelObj.ComputeTBDIS;
MRA.RunData.TBDIS  = MRA.ModelObj.TBDIS;
MRA.RunData.TBDISE = sqrt(MRA.ModelObj.TBDIS);
%%
MRA.chi2 = 'chi2CM'; MRA.ComputeCM;
MRA.fixPar = '1 5 6';
MRA.Fit;
