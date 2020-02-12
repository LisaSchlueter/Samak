MC = McRunGenerator('McFlag','Fake','InitFile',@ref_FakeRun_KNM1_StackingStudy,'nRuns',10,'qUfrac_RelErr',0.05,'qU_RelErr',0,'TimeSec_RelErr',0.2);
MC.ComputeFakeRun;
%%
%R = RunAnalysis('DataType','Fake','FakeStudyName','FakeRun_KNM1_StackingStudy','RunNr',1);
M = MultiRunAnalysis('RunList',1:10,'DataType','Fake','FakeStudyName','FakeRun_KNM1_StackingStudy','fixPar','5 6 7 8 9 10 11');
M.Fit;


%%

