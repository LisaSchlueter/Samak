%% 1. Feature: Create Basic Fake Run
BasicInitFile = @ref_BasicFakeRun_KNM1; 

MC = McRunGenerator('McFlag','Fake','InitFile',BasicInitFile);
MC.ComputeFakeRun;

A = RunAnalysis('RunNr',1,'DataType','Fake','FakeStudyName',extractAfter(func2str(BasicInitFile),'ref_'));
A.Fit;



%% 2. Feature: Create many runs with different properties: e.g. with qU variations

% settings
CleanUp = 'OFF';                                       % all previously calculated runs with same Init file are deleted
StudyInitFile = @ref_FakeRun_KNM1_StackingStudy_qUErr; % always create a new Init file for your specific study

% Get qU uncertainty from real data
Real = MultiRunAnalysis('RunList','KNM1','DataType','Real');
qU_RelErr = std(Real.SingleRunData.qU,0,2)./mean(Real.SingleRunData.qU,2);

% compute fake MC runs
nRuns = 10%Real.nRuns;
MC = McRunGenerator('McFlag','Fake','InitFile',StudyInitFile,'nRuns',nRuns,'qU_RelErr',max(qU_RelErr),'VarDist','Uniform');
if strcmp(CleanUp,'ON') % all previously calculated runs with same Init file are deleted
MC.CleanUpFakeRuns;
end
MC.ComputeFakeRun;

% analyze fake data
M = MultiRunAnalysis('RunList',1:nRuns,'DataType','Fake','FakeStudyName',extractAfter(func2str(StudyInitFile),'ref_'),'fixPar','5 6 7 8 9 10 11');
M.exclDataStart = 17;
M.Fit;