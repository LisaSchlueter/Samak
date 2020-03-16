% KNM2 Figure skating twins
range = 40;
E0 = knm2FS_GetE0Twins('SanityPlot','OFF');
RecomputeFlag = 'OFF';

%% load or calc
savedir = [getenv('SamakPath'),'knm2ana/knm2_FigureSkating/results'];
savename = sprintf('%sknm2FS_UniformFit_%.0feV.mat',savedir,range);

if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename)
else
    RunAnaArg = {'RunList','KNM2_Prompt',... % all KNM2 golden runs
        'fixPar','mNu E0 Bkg Norm',...           % free Parameter !!
        'DataType','Twin',...
        'FSDFlag','BlindingKNM2',...       % final state distribution (theoretical calculation)
        'ELossFlag','KatrinT2',...         % energy loss function     ( different parametrizations available)
        'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
        'chi2','chi2Stat',...              % statistics only
        'NonPoissonScaleFactor',1,...
        'TwinBias_Q',E0,...
        'ROIFlag','14keV'};
    
    %% build object of MultiRunAnalysis class
    A = MultiRunAnalysis(RunAnaArg{:});
    A.exclDataStart = A.GetexclDataStart(range);
    
    A.ModelObj.recomputeRF = 'ON';
    A.ModelObj.InitializeRF
    %% fit without corrections
    A.Fit;
    FitResult_ref  = A.FitResult;
    %% fit with broadening of RF + broadening/shift of FSD
    TimeSec = zeros(3,1);
    TimeSec(1) = sum(A.SingleRunData.TimeSec(1:171));
    TimeSec(2) = sum(A.SingleRunData.TimeSec(172:268));
    TimeSec(3) = sum(A.SingleRunData.TimeSec(269:361));
    MultiWeights = TimeSec./sum(TimeSec);
    MultiPos = [E0(1),E0(end-120),E0(end)]';
    MultiPosRel = MultiPos-wmean(MultiPos,MultiWeights);
    %%
    Sigma = std(E0).^2;
    A.ModelObj.LoadFSD('MultiPos',MultiPosRel,'MultiWeight',MultiWeights,...
        'SanityPlot','ON','Sigma',Sigma);
    A.ModelObj.ComputeTBDDS; A.ModelObj.ComputeTBDIS;
    MACE_Sigma = std(A.SingleRunData.qU,0,2);
    A.ModelObj.MACE_Sigma = MACE_Sigma;
    
    A.ModelObj.recomputeRF = 'ON';
    A.ModelObj.InitializeRF;
    A.Fit;
    FitResult_imp  = A.FitResult;
    save(savename,'FitResult_imp','FitResult_ref','E0','MACE_Sigma');
end
%% result
fprintf('---------------------------\n')
fprintf('mNuSq = %.3f eV^2  (ref) \n',FitResult_ref.par(1))
fprintf('mNuSq = %.3f eV^2  (imp) \n',FitResult_imp.par(1))
fprintf('---------------------------\n')


