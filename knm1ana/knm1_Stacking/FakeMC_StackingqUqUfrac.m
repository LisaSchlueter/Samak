% investigate effect of different qUfrac + TimeSec on stacking with fake MC
% test stacking qU correction
CleanUp = 'OFF';
exclDataStart = 14;
InitFile = @ref_FakeRun_KNM1_StackingStudy_qUqUfracErr;

%% Get qU uncertainty from real data
Real = MultiRunAnalysis('RunList','KNM1','DataType','Real');
qUfrac_RelErr = std(Real.SingleRunData.qUfrac,0,2)./mean(Real.SingleRunData.qUfrac,2);
qUfrac_RelErr = qUfrac_RelErr(2:end); % exclude -200eV point
qU_RelErr = std(Real.SingleRunData.qU,0,2)./mean(Real.SingleRunData.qU,2);
Time_RelErr   = std(Real.SingleRunData.TimeSec)./mean(Real.SingleRunData.TimeSec);
nRuns = Real.nRuns;
%plot(mean(Real.SingleRunData.qU,2)-18575,qUfrac_RelErr,'--x'); xlabel('qU - 18575 (eV)'); ylabel('qUfrac rel std'); PrettyFigureFormat; grid on;
%%

MC = McRunGenerator('McFlag','Fake','InitFile',InitFile,'nRuns',nRuns,'qUfrac_RelErr',max(qUfrac_RelErr),...
    'TimeSec_RelErr',Time_RelErr,'VarDist','Uniform','qU_RelErr',max(qU_RelErr));

if strcmp(CleanUp,'ON')
MC.CleanUpFakeRuns;
end

MC.ComputeFakeRun;

%%
M = MultiRunAnalysis('RunList',1:nRuns,'DataType','Fake','FakeStudyName',extractAfter(func2str(InitFile),'ref_'),'fixPar','5 6 7 8 9 10 11');
M.exclDataStart = exclDataStart;
M.Fit;
