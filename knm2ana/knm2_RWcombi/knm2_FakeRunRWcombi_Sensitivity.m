% Neutrino Mass Sensitivity from Fake MC runs
% with and without modification of FSDs

Mode = '3P';
savedir  = [getenv('SamakPath'),'knm2ana/knm2_RWcombi/results/'];
savename = [savedir,sprintf('knm2_FakeRunRWcombi_Sensitivity_%s.mat',Mode)];
RecomputeFlag = 'OFF';

if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename)
else
    
    switch Mode
        case 'Drift'
            maxRWdrift = 0.4;
            FSDArg = {'Dist','Rect','Sigma',maxRWdrift};
        case 'Gauss'
            sigma = 0.25;
            FSDArg = {'Dist','Gauss','Sigma',sigma};
        case '3P'
            %MultiPos     = [0.2042    0.7042   -0.9958];
            MultiPos     = [0.02042    0.07042   -0.09958];
            MultiWeights = [0.3929    0.3084    0.2987];
            FSDArg = {'Dist','Gauss','Sigma',0.1,'MultiPos',MultiPos,'MultiWeights',MultiWeights};
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
    %% Sensitivity when there is NO plasma fluctuation in data and model
    InitFileCombi =  @ref_FakeRun_KNM2_CD84_308x2hours; %84 column density, 308 runs a 2 hours
    M = RunAnalysis(CommonArg{:},'FakeInitFile',InitFileCombi,'TwinBias_Q',18573.70);
    M.InitModelObj_Norm_BKG;%('Recompute','ON');
    M.ModelObj.ComputeTBDDS;
    M.ModelObj.ComputeTBDIS;
    M.RunData.TBDIS = M.ModelObj.TBDIS;
    M.Fit;
    LmNuSq_ref = M.FitResult.err(1);
    %% Sensitivity when there is plasma fluctuation in data and model
    M.ModelObj.LoadFSD(FSDArg{:});
    M.InitModelObj_Norm_BKG;
    M.ModelObj.ComputeTBDDS;
    M.ModelObj.ComputeTBDIS;
    M.RunData.TBDIS = M.ModelObj.TBDIS;
    M.Fit;
    LmNuSq_plasma = M.FitResult.err(1);
    save(savename,'LmNuSq_ref','LmNuSq_plasma','FSDArg');
end
fprintf('-------------------------------------\n');
fprintf('Sensitivity nu-mass squared = %.2f meV^2 (no plasma time evolution) \n',LmNuSq_ref*1e3);
fprintf('Sensitivity nu-mass squared = %.2f meV^2 (with plasma time evolution) \n',LmNuSq_plasma*1e3);
fprintf('Relative increase =  %.1f %% \n',(LmNuSq_plasma-LmNuSq_ref)./LmNuSq_ref*100);






