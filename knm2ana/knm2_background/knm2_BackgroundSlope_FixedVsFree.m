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
    T.exclDataStop = 37; % exclude first background point
    T.chi2 = 'chi2CMShape';
    T.fixPar = 'mNu E0 Bkg Norm';
    T.InitFitPar; 
    FitParStd = 13.3.*1e-06; %std(FitPar(12,:),0,2); % standard deviation of background slope over samples
    
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
        T.RunData.TBDISE = sqrt(TBDIS(:,i));
        
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


