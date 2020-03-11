
E0 = knm2FS_GetE0Twins('SanityPlot','OFF');

RunAnaArg = {'RunList','KNM2_Prompt',... % all KNM2 golden runs
    'fixPar','E0 Bkg Norm',...           % free Parameter !!
    'DataType','Twin',...
    'FSDFlag','BlindingKNM2',...       % final state distribution (theoretical calculation)
    'ELossFlag','KatrinT2',...         % energy loss function     ( different parametrizations available)
    'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
    'chi2','chi2Stat',...              % statistics only
    'NonPoissonScaleFactor',1,...
    'TwinBias_Q',E0,...
    'ROIFlag','14keV'};

% build object of MultiRunAnalysis class
A = MultiRunAnalysis(RunAnaArg{:});

%%
% load runwise response functions
RFpath = [getenv('SamakPath'),'inputs/ResponseFunction/samakRF'];
RFrunwise = zeros(A.nRuns,A.ModelObj.nTe,A.ModelObj.nqU);

for i=1:A.nRuns
    RFfile =  sprintf('_Uniform_%.5gcm2_IsX%s_NIS%.0d_Bm%0.3gT_Bs%0.3gT_Ba%0.4gG_Temin%.3f_Temax%.3f_Bin%.0d_%s.mat',...
        A.SingleRunData.WGTS_CD_MolPerCm2(i),'Edep',A.NIS,A.SingleRunData.MACE_Bmax_T(i),A.SingleRunData.WGTS_B_T(i),...
        A.SingleRunData.MACE_Ba_T(i)*1e4,A.ModelObj.TeMin(i),A.ModelObj.TeMax,A.ModelObj.nTeBinningFactor,A.ModelObj.ELossFlag);
    
    RFrunwise(i,:,:) = load(RFfile,'RF');
    
end

