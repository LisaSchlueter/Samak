
A = RunAnalysis('FakeInitFile',@ref_FakeRun_KNM2_CD84_308x2hours,...
    'RunNr',1,'DataType','Fake',...
    'MosCorrFlag','OFF','ROIFlag','Default',...
    'fixPar','mNu E0 Bkg Norm');

A.ModelObj.recomputeRF = 'ON';

A.ModelObj.InitializeRF;%ComputeRF('ELossBinStep',0.04);
RFref = A.ModelObj.RF;
A.ModelObj.ComputeTBDDS;
A.ModelObj.ComputeTBDIS;
TBDIS = A.ModelObj.TBDIS;
%%
A.ModelObj.InitializeRF; % change eloss binning in computeRF by hand...
RFdiff = A.ModelObj.RF;
A.RunData.TBDIS = TBDIS;
A.Fit
