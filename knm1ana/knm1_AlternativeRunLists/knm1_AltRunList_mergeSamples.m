% merge several sub-files

savedir = [getenv('SamakPath'),'knm1ana/knm1_AlternativeRunLists/results/'];
nPixList = 1000;
FitResults = cell(nPixList,1);
RunLists   = cell(nPixList,1);
chi2 = 'chi2Stat';
NP = 1.064;

nSamples= [1000,1001];

savename1 = [savedir,sprintf('RandomHalfRunList_Unblinded_%s_NP%2g_%.0ffits.mat',chi2,NP,nSamples(1))];
savename2 = [savedir,sprintf('RandomHalfRunList_Unblinded_%s_NP%2g_%.0ffits.mat',chi2,NP,nSamples(2))];

d1 = importdata(savename1);
d2 = importdata(savename2);

FitResults = [d1.FitResults; d2.FitResults];
RunLists   = [d1.RunLists; d2.RunLists];
FitResults = {FitResults{1:sum(nSamples)-1}}'; % remove last sample
RunLists = {RunLists{1:sum(nSamples)-1}}';     % remove last sample
M = d1.M;
RunAnaArg = d1.RunAnaArg;

savenameNew = [savedir,sprintf('RandomHalfRunList_Unblinded_%s_NP%2g_%.0ffits_merged.mat',chi2,NP,numel(FitResults))];
save(savenameNew,'FitResults','M','RunAnaArg','RunLists','nSamples');








