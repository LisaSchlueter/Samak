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
%% plot neutrino mass w and w/o bkg slope
GetFigure;
h1 = histogram(FitPar(1,:)-FitParCM(1,:),'FaceColor',rgb('DodgerBlue'),'FaceAlpha',1);
xlabel(sprintf('{\\itm}_{{\\itB} slope free}^2 - {\\itm}_{{\\itB} slope fixed}^2 (eV^2)'))
ylabel('Occurence')
leg = legend(sprintf('Mean = %.1g eV^2',mean(FitPar(1,:)-FitParCM(1,:))),...
    'EdgeColor',rgb('Silver'),'Location','northwest');
PrettyFigureFormat('FontSize',22)
t = title(sprintf('{\\itB} slope constraint in CovMat = %.1f mcps / keV',FitParStd.*1e6),...
    'FontWeight','normal','FontSize',get(gca,'FontSize'));
plotname2 = sprintf('%sknm2_BkgSlope_FixedVsFree_mNuSq.png',plotdir);
print(plotname2,'-dpng','-r450');
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
