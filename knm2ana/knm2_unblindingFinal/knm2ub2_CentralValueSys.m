dir = [getenv('SamakPath'),'tritium-data/sensitivity/Knm2/'];
Blinding = 'OFF';
AnaFlag = 'Ring';%StackPixel';
Total = 1; % 0=no total, 1 = with total
%SysEffectsAll   = {'FSD','RF_BF','RF_RX','LongPlasma','FPDeff','NP','BkgPT','BkgPT','Bkg'}; %Bkg has to be last
SysEffectsAll   = {'FSD','RF_BF','RF_RX','LongPlasma','TASR','FPDeff','NP','BkgPT','Bkg'}; %Bkg has to be last

mNuSq    = zeros(numel(SysEffectsAll)+2+Total,1);
mNuSqErr = zeros(numel(SysEffectsAll)+2+Total,1);


if strcmp(AnaFlag,'Ring')
    nameStat = sprintf('%sRunSensAsimov_KNM2_Prompt_40eV_chi2Stat_MinuitMinosFit_Real_RingFull_qUOffsets.mat',dir);
else
    nameStat = sprintf('%sRunSensAsimov_KNM2_Prompt_39eV_chi2Stat_MinuitMinosFit_Real.mat',dir);
end
d = importdata(nameStat);
mNuSq(1)    = d.FitResult.par(1);
mNuSqErr(1) =  (-d.FitResult.errNeg(1)+d.FitResult.errPos(1))/2;
BptLoadLog = 0;

for i=1:numel(SysEffectsAll)
    if strcmp(AnaFlag,'Ring')
        nametmp = sprintf('%sRunSensAsimov_KNM2_Prompt_40eV_chi2CMShape_budget41_MinuitMinosFit_Real_%s_RingFull_qUOffsets.mat',dir,SysEffectsAll{i});
    else
        nametmp = sprintf('%sRunSensAsimov_KNM2_Prompt_39eV_chi2CMShape_budget40_MinuitMinosFit_Real_%s.mat',dir,SysEffectsAll{i});

    end
dtmp = importdata(nametmp);
mNuSq(i+1)    = dtmp.FitResult.par(1);
mNuSqErr(i+1) =  (-dtmp.FitResult.errNeg(1)+dtmp.FitResult.errPos(1))/2;

if strcmp(SysEffectsAll{i},'BkgPT') && BptLoadLog==0
    BptLoadLog=1;
elseif strcmp(SysEffectsAll{i},'BkgPT') && BptLoadLog==1
     namePullPT = sprintf('%sknm2ana/knm2_PngBkg/results/knm2ubfinal_Fit_Bpng-3.0mucpsPers_Real_40eV_mNuE0BkgNormBkgPTSlope_chi2Stat_StackPixel_KNM2_pull24.mat',getenv('SamakPath'));
     dtmp    = importdata(namePullPT);
     mNuSq(i+1) =dtmp.FitResult.par(1);
     mNuSqErr(i+1) =  (-dtmp.FitResult.errNeg(1)+dtmp.FitResult.errPos(1))/2;
end

end

%% pull term: fit with bkg slope as nuissance parameter with pull term
if ~strcmp(AnaFlag,'Ring')
    if any(ismember(SysEffectsAll,'Bkg'))
        namePull = sprintf('%sknm2ana/knm2_PngBkg/results/knm2ubfinal_FitBkgSlope_Real_40eV_mNuE0BkgNormBkgSlope_chi2Stat_StackPixel_BkgPull10_4.74mcpskeV.mat',getenv('SamakPath'));
        plotname = sprintf('%sknm2ana/knm2_PngBkg/plots/CentralValueSys_Real_StackedPixel_Stat.png',getenv('SamakPath'));
    else
        % namePull = sprintf('%sknm2ana/knm2_unblindingFinal/results/BestFit/knm2ub2_FitBkgSlope_Real_40eV_mNuE0BkgNormBkgSlope_chi2CMShape_StackPixel_BkgPull10_4.74mcpskeV_SysBudget38.mat',getenv('SamakPath'));
        namePull = sprintf('%sknm2ana/knm2_unblindingFinal/results/BestFit/knm2ub2_FitBkgSlope_Real_40eV_mNuE0BkgNormBkgSlope_chi2Stat_StackPixel_BkgPull10_4.74mcpskeV.mat',getenv('SamakPath'));
        plotname = sprintf('%sknm2ana/knm2_unblindingFinal/plots/Fit_BslopePull_Real_StackedPixel_Stat.png',getenv('SamakPath'));
    end
    dpull    = importdata(namePull);
    mNuSq(end-Total) = dpull.FitResult.par(1);
    mNuSqErr(end-Total) = (-dpull.FitResult.errNeg(1)+dpull.FitResult.errPos(1))/2;
    x = 1:numel(SysEffectsAll)+2+Total;    
else
    mNuSq(end-Total) = [];
    mNuSqErr(end-Total) = [];
    x = 1:numel(SysEffectsAll)+1+Total;
    plotname = sprintf('%sknm2ana/knm2_unblindingFinal/plots/Fit_BslopePull_Real_StackedPixel_Stat_MultiRing.png',getenv('SamakPath'));
end

if Total== 1
    if ~strcmp(AnaFlag,'Ring')
        nameTotal = sprintf('%sknm2ana/knm2_PngBkg/results/knm2ubfinal_Fit_Bpng-3.0mucpsPers_Real_40eV_mNuE0BkgNorm_chi2CMShape_StackPixel_KNM2_SysBudget40.mat',getenv('SamakPath'));
    else
        nameTotal = sprintf('%sknm2ana/knm2_PngBkg/results/knm2ubfinal_Fit_Bpng-3.0mucpsPers_Real_40eV_mNuE0BkgNormqU_chi2CMShape_Ring_KNM2_SysBudget41.mat',getenv('SamakPath'));
    end
    dtot = importdata(nameTotal);
    mNuSq(end)    =  dtot.FitResult.par(1);
    mNuSqErr(end) =  (-dtot.FitResult.errNeg(1)+dtot.FitResult.errPos(1))/2;
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
    ppull = plot(x(end-Total),mNuSq(end-Total)+rndNum,'.','MarkerSize',26,'LineWidth',1.5,'Color',rgb('Orange'));
    ppull2 = plot(x(end-2-Total),mNuSq(end-2-Total)+rndNum,'.','MarkerSize',26,'LineWidth',1.5,'Color',rgb('Orange'));
end
xticks(x); xlim([min(x)-0.5 max(x)+0.5]);
PrettyFigureFormat;
if ~strcmp(AnaFlag,'Ring')
    if any(ismember(SysEffectsAll,'BkgPT'))
        if Total== 1
            xticklabels({'Stat.','FSD','B-fields',sprintf('\\rhod\\sigma'),'Plasma',sprintf('\\epsilon_{FPD}'),sprintf('B_{NP}'),sprintf('B_{PT}'),sprintf('B_{PT}'),sprintf('B_{slope}'),sprintf('B_{slope}'),'Total'});
        else
            xticklabels({'Stat.','FSD','B-fields',sprintf('\\rhod\\sigma'),'Plasma',sprintf('\\epsilon_{FPD}'),sprintf('B_{NP}'),sprintf('B_{PT}'),sprintf('B_{PT}'),sprintf('B_{slope}'),sprintf('B_{slope}')});
        end
        legend([pstat,pcov,ppull],'stat. only','stat. and cov. mat.','stat. and pull term','Location','southwest','EdgeColor',rgb('Silver'));
    else
        xticklabels({'Stat.','FSD','B-fields',sprintf('\\rhod\\sigma'),'Plasma',sprintf('\\epsilon_{FPD}'),sprintf('B_{NP}'),sprintf('B_{slope}'),sprintf('B_{slope}')});
        legend([pstat,pcov,ppull],'stat. only','stat. and cov. mat.','stat. and pull term','Location','southwest','EdgeColor',rgb('Silver'));
    end
else
       if Total== 1
            xticklabels({'Stat.','FSD','B-fields',sprintf('\\rhod\\sigma'),'Plasma',sprintf('\\epsilon_{FPD}'),sprintf('B_{NP}'),sprintf('B_{PT}'),sprintf('B_{slope}'),'Total'});
       else
    xticklabels({'Stat.','FSD','B-fields',sprintf('\\rhod\\sigma'),'Plasma',sprintf('\\epsilon_{FPD}'),sprintf('B_{NP}'),sprintf('B_{PT}'),sprintf('B_{slope}')});
       end
    legend([pstat,pcov],'stat. only','stat. and cov. mat.','Location','southwest','EdgeColor',rgb('Silver'));
    
end
ylabel(sprintf('{\\itm}_\\nu^2 %s (eV^{ 2})',ystr));
ylim([min(mNuSq+rndNum)-0.01,max(mNuSq+rndNum)+0.01])
    %%
if strcmp(Blinding,'ON')
   plotname =  strrep(plotname,'.png','_Blinding.png');
end
if Total==1
    plotname =  strrep(plotname,'.png','_tot.png');
end
print(plotname,'-dpng','-r350');
fprintf('save plot to %s \n',plotname)