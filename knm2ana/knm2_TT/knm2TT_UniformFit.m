% KNM2 Figure skating twins
range = 40;
E0 = knm2FS_GetE0Twins('SanityPlot','OFF');
RecomputeFlag = 'OFF';
chi2 = 'chi2CMShape';
%% load or calc
savedir = [getenv('SamakPath'),'knm2ana/knm2_TT/results/'];
savename = sprintf('%sknm2TT_UniformFit_%.0feV_%s.mat',savedir,range,chi2);

if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename)
else
    fprintf('start fit\n');
    RunAnaArg = {'RunList','KNM2_Prompt',... % all KNM2 golden runs
        'fixPar','mNu E0 Bkg Norm',...       % free Parameter !!
        'DataType','Twin',...
        'FSDFlag','BlindingKNM2',...         % final state distribution (theoretical calculation) 
        'ELossFlag','KatrinT2',...           % energy loss function     ( different parametrizations available)
        'AnaFlag','StackPixel',...           % FPD segmentations -> pixel combination
        'chi2',chi2,...                      % statistics only
        'TwinBias_Q',E0,...
        'ROIFlag','14keV',...
        'SysBudget',34};
    if strcmp(chi2,'chi2Stat')
        RunAnaArg = [RunAnaArg,'NonPoissonScaleFactor',1];
    else
        RunAnaArg = [RunAnaArg,'NonPoissonScaleFactor',1.112];
    end
    %% build object of MultiRunAnalysis class
    A = MultiRunAnalysis(RunAnaArg{:});
    A.exclDataStart = A.GetexclDataStart(range);
    
    %% Broadeing with E0
    Sigma = std(E0);
    FSDArg = {'SanityPlot','ON','Sigma',Sigma};
    A.ModelObj.LoadFSD(FSDArg{:});
    A.ModelObj.ComputeTBDDS; A.ModelObj.ComputeTBDIS;

    A.Fit('CATS','ON');
    FitResult_imp  = A.FitResult;
    %save(savename,'FitResult_imp','FitResult_ref','E0','MACE_Sigma','A','FSDArg');
    save(savename,'FitResult_imp','FitResult_ref','E0','A','FSDArg');
end
