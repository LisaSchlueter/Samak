% KNM2 Figure skating twins
range = 40;
RecomputeFlag = 'OFF';

RFBinStep_i = 0.001;
RFBinStep = [0.001,0.005,0.01,0.02,0.04,0.06,0.08,0.1,0.2];

%% load or calc
savedir = [getenv('SamakPath'),'knm2ana/knm2_ResponseFunction/results/'];
savename = sprintf('%sFindOptStepSize_%.0fmeVRFi_min%.0fmeV_max%.0fmeV_%.0fSteps_%.0feV.mat',...
    savedir,1e3*RFBinStep_i,min(RFBinStep)*1e3,1e3*max(RFBinStep),numel(RFBinStep),range);

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
        'TwinBias_Q',18573.70,...
        'ROIFlag','14keV'};
    
    %% build object of MultiRunAnalysis class
    A = MultiRunAnalysis(RunAnaArg{:});
    A.exclDataStart = A.GetexclDataStart(range);
    
    A.ModelObj.RFBinStep = RFBinStep_i;
    A.ModelObj.InitializeRF;
    A.ModelObj.ComputeTBDDS;
    A.ModelObj.ComputeTBDIS;
    TBDIS_i = A.ModelObj.TBDIS;
    A.RunData.TBDIS = TBDIS_i;
    
    mNuSq = zeros(numel(RFBinStep),1);
    E0 = zeros(numel(RFBinStep),1);
    
    progressbar('find optimal rf step size')
    for i=1:numel(RFBinStep)
        progressbar(i/numel(RFBinStep));
        A.ModelObj.RFBinStep = RFBinStep(i);
        A.ModelObj.InitializeRF;
        A.Fit;
        mNuSq = A.FitResult.par(1);
        E0 = A.FitResult.par(2)+A.ModelObj.Q_i-18573.70; % endpoint offset from MC truth
    end
    
    MakeDir(savedir);
    save(savename,'mNuSq','E0','RFBinStep','RFBinStep_i','RunAnaArg','TBDIS_i')
end