

RunList = load('RunList100Good.mat');
RunList = RunList.RunList100Good;
RunList((RunList == 40610) | (RunList == 40612) | (RunList == 40539)) = [];


y = 1; n = 1;
for r = 1:length(RunList)
    
    R = load(num2str(RunList(r)));
    obj.StudyObject.TD = [num2str(RunList(r)),'ex2'];
    obj.StudyObject.WGTS_CD_MolPerCm2 = R.WGTS_CD_MolPerCm2;
    
    obj.WGTS_CD_MolPerCm2_RelErr = 0.08;
    obj.nTrials = 1000;
    obj.nRuns = 1;
    
    
    covmat_filename = sprintf('WGTSMACE_CovMat_%s_%g-molPercm2_%.2err_%gTrials_%u-Runs.mat',...
        obj.StudyObject.TD, obj.StudyObject.WGTS_CD_MolPerCm2,obj.WGTS_CD_MolPerCm2_RelErr,obj.nTrials,obj.nRuns);
    
    if exist(covmat_filename,'file')
        YesCovMat(y) = RunList(r);
        y = y + 1;
    else 
        NoCovMat(n) = RunList(r);
        n = n + 1;
    end
end