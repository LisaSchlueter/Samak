function CreateFakeRuns()
RunList = load('RunListFTFullCD.mat');
RunList = RunList.RunListFTFullCD;
RunList((RunList == 40610) | (RunList == 40612) | (RunList == 40539)) = [];

mpix = 'mpix'; % mpix
ringCutFlag = ''; %ex2

for r = 1:length(RunList)
    load([num2str(RunList(r)),mpix,ringCutFlag,'.mat']);
    A_MC = ref_RunAnalysis(RunList(r),mpix,ringCutFlag,'FPD_Segmentation','MULTIPIXEL',...
        'nTeBinningFactor',5,'FPD_Pixel',1:148);
    A_MC.ComputeTBDDS(); A_MC.ComputeTBDIS();
    A_MC.AddStatFluctTBDIS();
    
    TBDIS = A_MC.TBDIS;
    
    save(['fake-tritium-data/FAKE',num2str(RunList(r)),mpix,'.mat'],...
        'TBDIS','qU','TimeSec','qUfrac',...
        'WGTS_CD_MolPerCm2','WGTS_CD_MolPerCm2_SubRun',...
        'WGTS_MolFrac_TT','WGTS_MolFrac_DT','WGTS_MolFrac_HT',...
        'WGTS_MolFrac_TT_SubRun','WGTS_MolFrac_DT_SubRun','WGTS_MolFrac_HT_SubRun',...
        '-v7.3','-nocompression')
    
    % Save stacked pixels excluding two outer rings
    qU = mean(qU(:,1:124),2);
    TBDIS = sum(TBDIS(:,1:124),2);
    
    save(['fake-tritium-data/FAKE',num2str(RunList(r)),'ex2.mat'],...
        'TBDIS','qU','TimeSec','qUfrac',...
        'WGTS_CD_MolPerCm2','WGTS_CD_MolPerCm2_SubRun',...
        'WGTS_MolFrac_TT','WGTS_MolFrac_DT','WGTS_MolFrac_HT',...
        'WGTS_MolFrac_TT_SubRun','WGTS_MolFrac_DT_SubRun','WGTS_MolFrac_HT_SubRun',...
        '-v7.3','-nocompression')
    
end

end