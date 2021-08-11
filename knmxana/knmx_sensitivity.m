

FakeInitFile = @ref_FakeRun_TDRKatrin_1000days_10mcps;
A = RunAnalysis('RunNr',1,'DataType','Fake',...
    'fixPar','mNu E0 Norm Bkg',...
    'NonPoissonScaleFactor',1,...
    'FakeInitFile',FakeInitFile,...
    'ElossFlag','KatrinT2A20');

A.exclDataStart = A.GetexclDataStart(40);
A.Fit;
FitResultNom = A.FitResult; 
mNuSqStat90 = sqrt(FitResultNom.err(1).*1.64);

FakeInitFile0 = @ref_FakeRun_TDRKatrin_1000days_0mcps;
A0 = RunAnalysis('RunNr',1,'DataType','Fake',...
    'fixPar','mNu E0 Norm Bkg',...
    'NonPoissonScaleFactor',1,...
    'FakeInitFile',FakeInitFile0,...
    'ElossFlag','KatrinT2A20');

A0.exclDataStart = A0.GetexclDataStart(40);
A0.Fit;
FitResult_noB = A0.FitResult; 
mNuSqStat90_noB = sqrt(FitResult_noB.err(1).*1.64);
 

