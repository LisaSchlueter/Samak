DataType = 'Real';%Real';
AnaFlag = 'StackPixel';%Ring';
RingMerge = 'Full';%'None';

if strcmp(AnaFlag,'Ring')
    SysBudget = 39;
    if strcmp(RingMerge,'Full')
        AnaStr = AnaFlag;
    else
        AnaStr = sprintf('Ring%s',RingMerge);
    end
else
    SysBudget = 38;
    AnaStr = AnaFlag;
end
savedir = [getenv('SamakPath'),'knm2ana/knm2_background/results/'];
savename = sprintf('%sknm2_backgroundTimeEvolution140_%s_%s.mat',...
    savedir,DataType,AnaStr);

if exist(savename,'file') 
    load(savename,'FitResult','RunAnaArg','A');
else
    SigmaSq =  0.0124+0.0025;
    
    RunAnaArg = {'RunList','KNM2_Prompt',...
        'DataType',DataType,...
        'fixPar','mNu E0 Bkg Norm',...
        'RadiativeFlag','ON',...
        'minuitOpt','min ; minos',...
        'FSDFlag','BlindingKNM2',...
        'ELossFlag','KatrinT2A20',...
        'SysBudget',38,...
        'AnaFlag',AnaFlag,...
        'chi2','chi2Stat',...
        'TwinBias_Q',18573.7,...
        'NonPoissonScaleFactor',1,...
        'FSD_Sigma',sqrt(SigmaSq),...
        'TwinBias_FSDSigma',sqrt(SigmaSq),...
        'RingMerge',RingMerge};
    A = MultiRunAnalysis(RunAnaArg{:});
    A.exclDataStart = A.GetexclDataStart(40);
    %%
    Time = A.SingleRunData.qUfrac(end,:).*A.SingleRunData.TimeSec;
BkgPoint = A.SingleRunData.TBDIS(end,:)./Time;
LiveTime = hours(A.SingleRunData.StartTimeStamp-A.SingleRunData.StartTimeStamp(1));

   MakeDir(savedir);
   save(savename,'Time','BkgPoint','LiveTime');
end
%% Extract background for 140eV subrun
figure('Units','normalized','Position',[0.1,0.1,0.6,0.4]);
errorbar(LiveTime,BkgPoint,1.112.*sqrt(BkgPoint.*Time)./Time,...
    '.','MarkerSize',15,'LineWidth',2,'Color',rgb('Silver'),'MarkerEdgeColor',rgb('DodgerBlue'),'CapSize',0);
PrettyFigureFormat('FontSize',20);
xlabel('Time (hours)');
ylabel('Background rate (cps)');
%ylim([min(BkgPoint) max(BkgPoint)]);
xlim([-30 max(LiveTime)+30]);
title(sprintf('Backround subrun qU - 18574 = 140 eV'),'FontWeight','normal','FontSize',16);

[par, err, chi2min,dof,] = linFit(LiveTime',BkgPoint',1.112.*(sqrt(BkgPoint.*Time)./Time)');
hold on;
pFit = plot(LiveTime,p(1).*LiveTime+p(2),'-','LineWidth',2);
leg = legend(pFit,sprintf('Linear fit (stat + NP): %.2g \\pm %.1g mcps/h',par(1)*1e3,err(1)*1e3),...
    'FontSize',16,'EdgeColor',rgb('Silver'),'Location','northwest');
plotdir = strrep(savedir,'results','plots');
MakeDir(plotdir);
plotname = strrep(strrep(savename,'results','plots'),'.mat','.png');
print(plotname,'-dpng','-r350');
