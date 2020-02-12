% Neutrino Mass Sensitivity from Fake MC runs
% Investiation of plasma time evolution systematics
% -> modification of FSDs
% -> covariance matrix
Mode = 'Gauss';
SysErr = 0.04;
savedir  = [getenv('SamakPath'),'knm2ana/knm2_RWcombi/results/'];
savename = [savedir,sprintf('knm2_FakeRunRWcombi_Systematics_%s_%.0fmeVErr.mat',Mode,SysErr*1e3)];
RecomputeFlag = 'OFF';

if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename,'LmNuSq_cm','LmNuSq_stat','FitResult_stat','FitResult_cm','FSDArg','CM')
else
    
    switch Mode
        case 'Drift'
            Sigma = 0.4; % max RW drift
            Dist = 'Rect';
        case 'Gauss'
             SigmaErr = SysErr;
             MultiPosErr = 0;
            Sigma = 0.25;
            Dist='Gauss';
            FSDArg = {'Sigma',Sigma,'Dist',Dist};
            MultiPos = '';
            MultiWeights = '';
        case '3P'
            MultiPosErr = SysErr;
            MultiPos     = [0.02042    0.07042   -0.09958];
            MultiWeights = [0.3929    0.3084    0.2987];            
            Sigma = 0.001;
            Dist='Gauss';
            FSDArg = {'Sigma',Sigma,'MultiPos',MultiPos,'MultiWeights',MultiWeights,'Dist',Dist};
            %CMArg = {'MultiPos',MultiPos,'MultiWeights',MultiWeights,'Sigma',Sigma};
            SigmaErr = 0;
    end
    
    CommonArg = {'RunNr',1,...% has no meaning
        'DataType','Fake',...
        'FSDFlag','BlindingKNM2',...
        'ELossFlag','KatrinT2',...
        'exclDataStart',11,... % 11==40eV range (28 subruns)
        'chi2','chi2Stat',...
        'RingMerge','Full',...
        'minuitOpt','min;minos',...
        'NonPoissonScaleFactor',1,...
        'AnaFlag','StackPixel',...
        'fixPar','mNu E0 Bkg Norm',...
        'RingList',1:12};
    
    %% Init Model and MC data with plasma evolution
    InitFileCombi =  @ref_FakeRun_KNM2_CD84_308x2hours; %84 column density, 308 runs a 2 hours
    M = RunAnalysis(CommonArg{:},'FakeInitFile',InitFileCombi,'TwinBias_Q',18573.70);
    
    M.ModelObj.FSD_Sigma = Sigma;
    M.ModelObj.FSD_MultiPos = MultiPos;
    M.ModelObj.FSD_MultiWeights = MultiWeights;
    M.ModelObj.FSD_Dist = Dist;
    M.ModelObj.LoadFSD;
    M.InitModelObj_Norm_BKG;
    M.ModelObj.ComputeTBDDS;
    M.ModelObj.ComputeTBDIS;
    M.RunData.TBDIS = M.ModelObj.TBDIS;
    
    %% Stat only
    M.chi2 = 'chi2Stat';
    M.Fit;
    FitResult_stat  = M.FitResult;
    LmNuSq_stat     = M.FitResult.err(1);
    
    %% Stat + syst
    M.chi2 = 'chi2CMShape';


    M.ComputeCM('SysEffects',struct('RW','ON'),...
        'BkgCM','OFF',...
        'nTrials',1000,...
        'RW_MultiPosErr',MultiPosErr,...
        'RW_SigmaErr',SigmaErr);
    M.Fit;
    FitResult_cm  = M.FitResult;
    LmNuSq_cm     = M.FitResult.err(1);
    
   %% save
    CM = M.FitCM_Obj; %covariance matrix object 
    save(savename,'LmNuSq_cm','LmNuSq_stat','FitResult_stat','FitResult_cm','FSDArg','CM');
end
fprintf('Sys Err = %.0f meV -------------------------------------\n',SysErr*1e3);
fprintf('Stat only:     %.2f meV^2  \n',LmNuSq_stat*1e3);
fprintf('Stat. + syst:  %.2f meV^2  \n',LmNuSq_cm*1e3);
fprintf('Syst only:  %.2f meV^2  \n',sqrt(LmNuSq_cm^2-LmNuSq_stat^2)*1e3);
%% plot covariance matrix and correlation matrix
%CM.PlotCM('qUWindowIndexMax',10,'saveplot','ON','savename',sprintf('_%.0fmeVErr',MultiPosErr*1e3));
% CM.PlotCorr('qUWindowIndexMax',10,'saveplot','ON');

