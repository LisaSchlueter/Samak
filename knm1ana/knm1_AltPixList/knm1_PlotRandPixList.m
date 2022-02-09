% Choose random pixellist with half of KNM1 pixels
savedir = [getenv('SamakPath'),'knm1ana/knm1_AltPixList/results/'];
nPixList = 1000;
FitResults = cell(nPixList,1);
PixLists   = cell(nPixList,1);
chi2 = 'chi2Stat';
NP = 1.064;

savename = [savedir,sprintf('RandomHalfPixList_Unblinded_%s_NP%2g_%.0ffits.mat',chi2,NP,nPixList)];
%%
if exist(savename,'file')
    d = importdata(savename);
     fprintf('load %s \n',savename);
else
    fprintf(2,'file not found %s \n',savename);
end

mNuSq = cellfun(@(x) x.par(1),d.FitResults);
histogram(mNuSq);