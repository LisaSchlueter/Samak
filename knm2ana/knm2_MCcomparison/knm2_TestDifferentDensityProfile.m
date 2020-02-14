

saveStrResult =  [getenv('SamakPath'),'/knm2ana/knm2_MCcomparison/results/ISProb_DifferentCDprofiles.mat'];

if exist(saveStrResult,'file') 
    load(saveStrResult);
else
    T = ref_FakeRun_KNM2_RFcomparison('reComputeRF','OFF');
    ISProb_Table = T.ComputeISProb('WGTS_DensityProfile','Table','Method','New');
    ISProb_File  = T.ComputeISProb('WGTS_DensityProfile','File','Method','New');
    ISProb_Flat  = T.ComputeISProb('WGTS_DensityProfile','Flat','Method','New');
    save(saveStrResult,'ISProb_File','ISProb_Flat','ISProb_Table');
end