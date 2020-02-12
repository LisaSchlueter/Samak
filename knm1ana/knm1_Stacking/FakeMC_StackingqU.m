% investigate effect of different qU on stacking with fake MC
% test stacking qU correction

%% settings
CleanUp = 'OFF';
exclDataStart = 14;
InitFile = @ref_FakeRun_KNM1_StackingStudy_qUErr;

%% Get qU uncertainty from real data
Real = MultiRunAnalysis('RunList','KNM1','DataType','Real');
qU_RelErr = std(Real.SingleRunData.qU,0,2)./mean(Real.SingleRunData.qU,2);
%plot(mean(Real.SingleRunData.qU,2)-18575,qU_RelErr,'--x'); xlabel('qU - 18575 (eV)'); ylabel('qU rel std'); PrettyFigureFormat; grid on;
%% compute fake MC runs
nRuns = Real.nRuns;
MC = McRunGenerator('McFlag','Fake','InitFile',InitFile,'nRuns',nRuns,'qU_RelErr',max(qU_RelErr),'VarDist','Uniform');

if strcmp(CleanUp,'ON')
MC.CleanUpFakeRuns;
end

MC.ComputeFakeRun;

%%
M = MultiRunAnalysis('RunList',1:nRuns,'DataType','Fake','FakeStudyName',extractAfter(func2str(InitFile),'ref_'),'fixPar','5 6 7 8 9 10 11');
M.exclDataStart = exclDataStart;
M.Fit;

% %%
%  Mcorr = MultiRunAnalysis('RunList',1:nRuns,'DataType','Fake','FakeStudyName',extractAfter(func2str(InitFile),'ref_'),'fixPar','5 6 7 8 9 10 11',...
%                          'StackqUCorrFlag','ON','StackTolerance',10);
%  Mcorr.Fit;
