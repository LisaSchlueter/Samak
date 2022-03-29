% study to look at background slope
% (fixed to 0 + uncertainty) vs. free
% fit randomized MC twins

%% set up model
range = 40;
NonPoissonScaleFactor = 1.112;
nSamples = 1e3;

RunAnaArg = {'RunList','KNM2_Prompt',...
    'fixPar','mNu E0 Bkg Norm',...
    'DataType','Twin',...
    'FSDFlag','BlindingKNM2',...
    'ELossFlag','KatrinT2A20',...
    'AnaFlag','StackPixel',...
    'chi2','chi2Stat',...
    'ROIFlag','Default',...
    'SynchrotronFlag','ON',...
    'AngularTFFlag','ON',...
    'ISCSFlag','Edep',...
    'TwinBias_Q',18573.56,...
    'NonPoissonScaleFactor',NonPoissonScaleFactor};

savedir = [getenv('SamakPath'),'knm2ana/knm2_background/results/'];
MakeDir(savedir);

savenameFree = sprintf('%sknm2_BackgroundSlope_Free_%.0fSamples.mat',savedir,nSamples);
if exist(savenameFree,'file')
    load(savenameFree)
else
    T = MultiRunAnalysis(RunAnaArg{:});
    T.exclDataStart = T.GetexclDataStart(range);
    
    %% fit randomized MC with background slope
    T.fixPar = 'mNu E0 Bkg Norm BkgSlope';
    T.InitFitPar;
    [FitPar, FitErr, FitChi2min, dof, TBDIS_mc]  = T.FitTwin('nSamples',nSamples);
    save(savenameFree,'FitPar','FitErr','FitChi2min','dof','TBDIS_mc','RunAnaArg');
end

%% fit ranomized MC with background slope systematics
savenameFix = sprintf('%sknm2_BackgroundSlope_FixCM_%.0fSamples.mat',savedir,nSamples);

if exist(savenameFix,'file')
    load(savenameFix)
else
    T = MultiRunAnalysis(RunAnaArg{:});
    T.exclDataStart = T.GetexclDataStart(range);
  
    T.chi2 = 'chi2CMShape';
    T.fixPar = 'mNu E0 Bkg Norm';
    T.InitFitPar; 
    FitParStd =  std(FitPar(12,:),0,2); % standard deviation of background slope over samples
    
    CmArg = {'BkgCM','ON',...
        'SysEffects',struct('BkgShape','ON'),...
        'MaxSlopeCpsPereV',FitParStd,...
        'BkgMod','Gauss'};
    T.ComputeCM(CmArg{:});
    
    % init
    FitParCM     = zeros(T.nPar,nSamples);
    FitErrCM     = zeros(T.nPar,nSamples);
    FitChi2minCM = zeros(nSamples,1);
    dofCM        = 0;
    
    
    for i=1:nSamples
        progressbar(i/nSamples)
        
        T.RunData.TBDIS = TBDIS_mc(:,i); % use the sampe samples!
        T.RunData.TBDISE = sqrt(TBDIS_mc(:,i));
        
        T.Fit;
        FitParCM(:,i)   = T.FitResult.par;
        FitErrCM(:,i)   = T.FitResult.err;
        FitChi2minCM(i) = T.FitResult.chi2min;
        dofCM           = T.FitResult.dof;
        
        T.ModelObj.SetFitBias(0);
    end
    save(savenameFix,'FitParCM','FitErrCM','FitChi2minCM','dofCM','TBDIS_mc','RunAnaArg','CmArg','FitParStd',...
                     'FitPar','FitErr','FitChi2min','dof');
end

%% plots result with free bkg slope
plotdir = strrep(savedir,'results','plots');
MakeDir(plotdir);
plotname1 = sprintf('%sknm2_BkgSlope_FixedVsFree_ScatterHist.png',plotdir);
ScatterHist2(FitPar(1,:),FitPar(12,:).*1e6,'RefLine','ON',...
    'xName',sprintf('{\\itm}_\\nu^2 (eV^2)'),'yname',sprintf('{\\it B} slope (mcps / keV)'),...
    'SaveAs',plotname1);


%% scatter plot m2, bkg slope
CorrMat = corrcoef(FitPar(1,:),FitPar(12,:));
CorrMat_i = corrcoef(FitPar'); % correlation matrix

nReSample = 5;%1e4;

RandIdx = randi(size(FitPar,2),[size(FitPar,2),nReSample]);
CorrMat_sample = zeros(2,2,nReSample);

for i=1:nReSample
    CorrMat_sample(:,:,i) = corrcoef(FitPar([1,12],RandIdx(:,i)')');
end
MeanCorrMat = mean(CorrMat_sample,3);
StdCorrMat = std(CorrMat_sample,0,3);

%%

f1 = figure('Units','normalized','Position',[0.1,0.1,0.4,0.4]);
LineArg = {':','LineWidth',2,'Color',rgb('Gray')};
plot(linspace(-5,5,1e2),linspace(1e6.*median(FitPar(12,:)),1e6.*median(FitPar(12,:)),1e2),LineArg{:});
hold on;
plot(linspace(median(FitPar(1,:)),median(FitPar(1,:)),1e2),linspace(-50,50,1e2),LineArg{:});
p = dscatter(FitPar(1,:)',FitPar(12,:)'.*1e6);
pnone = plot(NaN,NaN,'wo','MarkerFaceColor','none','MarkerEdgeColor','none');

leg = legend(pnone,sprintf('\\rho = %.2f \\pm %.2f',CorrMat(2),StdCorrMat(2)),'Location','northwest');
PrettyLegendFormat(leg);
leg.ItemTokenSize = [0,0];
leg.TextColor = rgb('DeepPink');
leg.FontSize = 16;

xlabel(sprintf('{\\itm}_\\nu^{ 2} (eV^{ 2})'));
ylabel(sprintf('{\\its}_{qU} (mcps/keV)'));
PrettyFigureFormat('FontSize',18);

%t = text(-1.46,3,'median','FontSize',get(gca,'FontSize')+2,'Color',rgb('Gray'));
%t2 = text(0.06,48,'median','FontSize',get(gca,'FontSize')+2,'Rotation',270,'Color',rgb('Gray'));
xlim([-1.5 1.5]);
ylim([-50 50]);

plotname = sprintf('%sknm2_BkgSlope_FixedVsFree.pdf',plotdir);
export_fig(plotname);

%% plot neutrino mass w and w/o bkg slope
GetFigure;
h1 = histogram(FitPar(1,:)-FitParCM(1,:),'FaceColor',rgb('DodgerBlue'),'FaceAlpha',0.8,'EdgeColor',rgb('Blue'));
xlabel(sprintf('\\Delta{\\itm}_\\nu^2 (eV^{ 2})'));%sprintf('{\\itm}_\\nu^2 - {\\itm}_{{\\itB} slope fixed}^2 (eV^2)'))
ylabel('Occurence')
leg = legend(sprintf('Median = %.1g eV^2, std = %.1g eV^2',median(FitPar(1,:)-FitParCM(1,:)),std(FitPar(1,:)-FitParCM(1,:))),...
   'Location','northwest');
PrettyLegendFormat(leg);
PrettyFigureFormat('FontSize',22)
% t = title(sprintf('{\\itB} slope constraint in CovMat = %.1f mcps / keV',FitParStd.*1e6),...
%     'FontWeight','normal','FontSize',get(gca,'FontSize'));
plotname2 = sprintf('%sknm2_BkgSlope_FixedVsFree_mNuSq.pdf',plotdir);
export_fig(plotname2);
%print(plotname2,'-dpng','-r450');
%% plot uncertainties on neutrino mass
GetFigure;
ErrBfree  = FitErr(1,FitErr(1,:)>0.05 & FitErr(1,:)<0.8);
ErrBfixed = FitErrCM(1,FitErrCM(1,:)>0.05 & FitErrCM(1,:)<0.8);
ErrDiff = FitErr(1,:)-FitErrCM(1,:);
h1 = histogram(ErrBfree,'FaceColor',rgb('DodgerBlue'),'FaceAlpha',0.3);
hold on;
h2 = histogram(ErrBfixed,'FaceColor',rgb('Orange'),'FaceAlpha',0.3);
hold off;
xlabel(sprintf('{\\itm}_{err}^2 (eV^2)'))
ylabel('Occurence')
leg = legend([h1,h2],sprintf('{\\itB} slope free:   mean = %.2g eV^2',mean(ErrBfree)),...
    sprintf('{\\itB} slope fixed: mean = %.2g eV^2',mean(ErrBfixed)),...
    'EdgeColor',rgb('Silver'),'Location','northwest');
PrettyFigureFormat('FontSize',22)
xlim([0 0.7])
plotname3 = sprintf('%sknm2_BkgSlope_FixedVsFree_mNuSqErr.png',plotdir);
print(plotname3,'-dpng','-r450');

%%
%std(FitPar(1,:))
%std(FitParCM(1,:))
