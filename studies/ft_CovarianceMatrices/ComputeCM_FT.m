% Scrip to Compute CovarianceMatrices for all Runs
% To be computed on server
addpath(genpath('../../../Samak2.0'));
M = MultiRunAnalysis('RunList','StackCD100_3hours','chi2','chi2Stat','DataEffCor','RunSummary');
M = MultiRunAnalysis('RunList','StackCD100all','chi2','chi2Stat','DataEffCor','RunSummary');
RunList = M.SingleRunData.Runs(M.SingleRunData.Select_all);
for i=1:numel(RunList)
    fprintf(2,'RunNr %.0f',i)
    R = RunAnalysis('RunNr',RunList(i),'chi2','chi2Stat','DataEffCor','RunSummary','ELossFlag','Abdurashitov');
    R.ComputeCM('nTrials',1000);
end