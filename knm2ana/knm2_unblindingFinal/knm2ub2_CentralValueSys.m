dir = [getenv('SamakPath'),'tritium-data/sensitivity/Knm2/'];
Blinding = 'OFF';
AnaFlag = 'StackPixel';
SysEffectsAll   = {'FSD','RF_BF','RF_RX','LongPlasma','FPDeff','NP','Bkg'}; %Bkg has to be last
mNuSq    = zeros(numel(SysEffectsAll)+2,1);
mNuSqErr = zeros(numel(SysEffectsAll)+2,1);

if strcmp(AnaFlag,'Ring')
    nameStat = sprintf('%sRunSensAsimov_KNM2_Prompt_40eV_chi2Stat_MinuitMinosFit_Real_RingFull.mat',dir);
else
    nameStat = sprintf('%sRunSensAsimov_KNM2_Prompt_39eV_chi2Stat_MinuitMinosFit_Real.mat',dir);
end
d = importdata(nameStat);
mNuSq(1)    = d.FitResult.par(1);
mNuSqErr(1) =  (-d.FitResult.errNeg(1)+d.FitResult.errPos(1))/2;

for i=1:numel(SysEffectsAll)
    if strcmp(AnaFlag,'Ring')
nametmp = sprintf('%sRunSensAsimov_KNM2_Prompt_40eV_chi2CMShape_budget39_MinuitMinosFit_Real_%s_RingFull.mat',dir,SysEffectsAll{i});
    else
        nametmp = sprintf('%sRunSensAsimov_KNM2_Prompt_39eV_chi2CMShape_budget38_MinuitMinosFit_Real_%s.mat',dir,SysEffectsAll{i});

    end
dtmp = importdata(nametmp);
mNuSq(i+1)    = dtmp.FitResult.par(1);
mNuSqErr(i+1) =  (-dtmp.FitResult.errNeg(1)+dtmp.FitResult.errPos(1))/2;
end

%% pull term: fit with bkg slope as nuissance parameter with pull term
if ~strcmp(AnaFlag,'Ring')
   % namePull = sprintf('%sknm2ana/knm2_unblindingFinal/results/BestFit/knm2ub2_FitBkgSlope_Real_40eV_mNuE0BkgNormBkgSlope_chi2CMShape_StackPixel_BkgPull10_4.74mcpskeV_SysBudget38.mat',getenv('SamakPath'));
    namePull = sprintf('%sknm2ana/knm2_unblindingFinal/results/BestFit/knm2ub2_FitBkgSlope_Real_40eV_mNuE0BkgNormBkgSlope_chi2Stat_StackPixel_BkgPull10_4.74mcpskeV.mat',getenv('SamakPath'));
  
    dpull    = importdata(namePull);
    mNuSq(end) = dpull.FitResult.par(1);
    mNuSqErr(end) = (-dpull.FitResult.errNeg(1)+dpull.FitResult.errPos(1))/2;
    x = 1:numel(SysEffectsAll)+2;
    plotname = sprintf('%sknm2ana/knm2_unblindingFinal/plots/Fit_BslopePull_Real_StackedPixel_Stat.png',getenv('SamakPath'));
else
    mNuSq(end) = [];
    mNuSqErr(end) = [];
    x = 1:numel(SysEffectsAll)+1;
    plotname = sprintf('%sknm2ana/knm2_unblindingFinal/plots/Fit_BslopePull_Real_StackedPixel_Stat_MultiRing.png',getenv('SamakPath'));
end

switch Blinding
    case 'ON'
        rndNum = randn(1);
        ystr = '+ random number';
    case 'OFF'
        rndNum = 0;
        ystr = '';
end

f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.45]);

pcov = plot(x,mNuSq+rndNum,'.:','MarkerSize',26,'LineWidth',1.5,'Color',rgb('DodgerBlue'));
hold on;
pstat = plot(x(1),mNuSq(1)+rndNum,'k.','MarkerSize',26,'LineWidth',1.5);
if ~strcmp(AnaFlag,'Ring')
ppull = plot(x(end),mNuSq(end)+rndNum,'.','MarkerSize',26,'LineWidth',1.5,'Color',rgb('Orange'));
end
xticks(x); xlim([min(x)-0.5 max(x)+0.5]);
PrettyFigureFormat;
if ~strcmp(AnaFlag,'Ring')
xticklabels({'Stat.','FSD','B-fields',sprintf('\\rhod\\sigma'),'Plasma',sprintf('\\epsilon_{FPD}'),sprintf('B_{NP}'),sprintf('B_{slope}'),sprintf('B_{slope}')});
legend([pstat,pcov,ppull],'stat. only','stat. and cov. mat.','stat. and pull term','Location','southwest','EdgeColor',rgb('Silver'));
else
    xticklabels({'Stat.','FSD','B-fields',sprintf('\\rhod\\sigma'),'Plasma',sprintf('\\epsilon_{FPD}'),sprintf('B_{NP}'),sprintf('B_{slope}')});
legend([pstat,pcov],'stat. only','stat. and cov. mat.','Location','southwest','EdgeColor',rgb('Silver'));

end
ylabel(sprintf('{\\itm}_\\nu^2 %s (eV^{ 2})',ystr));
ylim([min(mNuSq(end)+rndNum)-0.01,max(mNuSq(1)+rndNum)+0.01])
%%
if strcmp(Blinding,'ON')
   plotname =  strrep(plotname,'.png','_Blinding.png');
end
print(plotname,'-dpng','-r350');
fprintf('save plot to %s \n',plotname)