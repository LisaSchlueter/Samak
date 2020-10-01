dir = [getenv('SamakPath'),'tritium-data/sensitivity/Knm2/'];

SysEffectsAll   = {'FSD','RF_BF','RF_RX','LongPlasma','FPDeff','NP','Bkg'}; %Bkg has to be last
mNuSq    = zeros(numel(SysEffectsAll)+2,1);
mNuSqErr = zeros(numel(SysEffectsAll)+2,1);
 
nameStat = sprintf('%sRunSensAsimov_KNM2_Prompt_39eV_chi2Stat_MinuitMinosFit_Real.mat',dir);
d = importdata(nameStat);
mNuSq(1)    = d.FitResult.par(1);
mNuSqErr(1) =  (-d.FitResult.errNeg(1)+d.FitResult.errPos(1))/2;

for i=1:numel(SysEffectsAll)
nametmp = sprintf('%sRunSensAsimov_KNM2_Prompt_39eV_chi2CMShape_budget38_MinuitMinosFit_Real_%s.mat',dir,SysEffectsAll{i});
dtmp = importdata(nametmp);
mNuSq(i+1)    = dtmp.FitResult.par(1);
mNuSqErr(i+1) =  (-dtmp.FitResult.errNeg(1)+dtmp.FitResult.errPos(1))/2;
end

%% pull term: fit with bkg slope as nuissance parameter with pull term
namePull = sprintf('%sknm2ana/knm2_unblinding1/results/Fit_BslopePull_Real_StackedPixel_Stat.mat',getenv('SamakPath'));
dpull    = importdata(namePull);
mNuSq(end) = dpull.FitResult.par(1);
mNuSqErr(end) = (-dpull.FitResult.errNeg(1)+dpull.FitResult.errPos(1))/2;

rndNum = randn(1);
f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.45]);
x = 1:numel(SysEffectsAll)+2;
pcov = plot(x,mNuSq+rndNum,'.:','MarkerSize',26,'LineWidth',1.5,'Color',rgb('DodgerBlue'));
hold on;
pstat = plot(x(1),mNuSq(1)+rndNum,'k.','MarkerSize',26,'LineWidth',1.5);
ppull = plot(x(end),mNuSq(end)+rndNum,'.','MarkerSize',26,'LineWidth',1.5,'Color',rgb('Orange'));
xticks(x); xlim([min(x)-0.5 max(x)+0.5]);
PrettyFigureFormat;
xticklabels({'Stat.','FSD','B-fields',sprintf('\\rhod\\sigma'),'Plasma',sprintf('\\epsilon_{FPD}'),sprintf('B_{NP}'),sprintf('B_{slope}'),sprintf('B_{slope}')});
ylabel(sprintf('{\\itm}_\\nu^2 + random number (eV^{ 2})'));
ylim([min(mNuSq(end)+rndNum)-0.01,max(mNuSq(1)+rndNum)+0.01])
legend([pstat,pcov,ppull],'stat. only','stat. and cov. mat.','stat. and pull term','Location','southwest','EdgeColor',rgb('Silver'));
%%
plotname = sprintf('%sknm2ana/knm2_unblinding1/plots/Fit_BslopePull_Real_StackedPixel_Stat.png',getenv('SamakPath'));
print(plotname,'-dpng','-r350');
