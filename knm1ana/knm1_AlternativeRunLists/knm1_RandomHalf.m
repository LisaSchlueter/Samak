% Choose random runlist with half of KNM1 runs
savedir = [getenv('SamakPath'),'knm1ana/knm1_AlternativeRunLists/results/'];
nRunList = 100;
FitResults = cell(nRunList,1);
RunLists   = cell(nRunList,1);
chi2 = 'chi2Stat';
savename = [savedir,sprintf('RandomHalfRunList_Unblinded_%s_%.0f.mat',chi2,nRunList)];

if exist(savename,'file')
    load(savename,'FitResults','RunLists');
else
progressbar('Random RunLists ');
for i=1:nRunList
    progressbar(i/nRunList);
M = MultiRunAnalysis('RunList','KNM1_Random','exclDataStart',14,'DataType','Real','chi2',chi2,...
    'fixPar','5 6 7 8 9 10 11');
M.Fit;

FitResults{i} = M.FitResult;
RunLists{i} = M.RunList;
end

save(savename,'FitResults','RunLists');
end
%% histogram neutrino mass and 1 sigma
Both = 'OFF';
if strcmp(Both,'ON')
    savedir = [getenv('SamakPath'),'knm1ana/knm1_AlternativeRunLists/results/'];
    
savename20 = [savedir,sprintf('RandomHalfRunList_%s_%.0f.mat',chi2,20)]; 
d20 = importdata(savename20);
mNuSq    = cell2mat(cellfun(@(x) x.par(1),d20.FitResults,'UniformOutput',0));
mNuSqErr = cell2mat(cellfun(@(x) x.err(1),d20.FitResults,'UniformOutput',0));

savename44 = [savedir,sprintf('RandomHalfRunList_%s_%.0f.mat',chi2,44)];
d44 = importdata(savename44);
mNuSq    = [mNuSq;cell2mat(cellfun(@(x) x.par(1),d44.FitResults,'UniformOutput',0))];
mNuSqErr = [mNuSqErr;cell2mat(cellfun(@(x) x.err(1),d44.FitResults,'UniformOutput',0))];


savename10 = [savedir,sprintf('RandomHalfRunList_%s_%.0f.mat',chi2,10)];
d10 = importdata(savename10);
mNuSq    = [mNuSq;cell2mat(cellfun(@(x) x.par(1),d10.FitResults,'UniformOutput',0))];
mNuSqErr = [mNuSqErr;cell2mat(cellfun(@(x) x.err(1),d10.FitResults,'UniformOutput',0))];

else
mNuSq    = cell2mat(cellfun(@(x) x.par(1),FitResults,'UniformOutput',0));
mNuSqErr = cell2mat(cellfun(@(x) x.err(1),FitResults,'UniformOutput',0));
end
plotdir = strrep(savedir,'results','plots');
plotname1 = [plotdir,sprintf('RandomHalfRunList_%s_%.0f_mNuSq.png',chi2,numel(mNuSq))];
plotname2 = [plotdir,sprintf('RandomHalfRunList_%s_%.0f_mNuSqErr.png',chi2,numel(mNuSq))];

f1 = figure(1);
set(f1, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.6, 0.6]);
h1 =histogram(mNuSq);
h1.FaceColor = rgb('SteelBlue'); h1.FaceAlpha = 0.8;
xlabel(sprintf('m_\\nu^2 (eV^2)'));
ylabel('# random runlist');
PrettyFigureFormat;
set(gca,'FontSize',22);
leg = legend(sprintf('%.0f samples',numel(mNuSq))); legend boxoff
print(f1,plotname1,'-dpng','-r450');


f2 = figure(2);
set(f2, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.6, 0.6]);
h2 =histogram(mNuSqErr);
h2.FaceColor = rgb('SteelBlue'); h2.FaceAlpha = 0.8;
xlabel(sprintf('\\sigma(m_\\nu^2) (eV^2)'));
ylabel('# random runlist');
PrettyFigureFormat;
set(gca,'FontSize',22);
leg = legend(sprintf('%.0f samples',numel(mNuSq))); legend boxoff
print(f2,plotname2,'-dpng','-r450');


