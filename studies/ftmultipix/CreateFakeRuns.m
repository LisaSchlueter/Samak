clear;

%D = load('40263.mat');

DX = ref_runsummaries(40263);
DX.TimeSec = DX.TimeSec;


for ii = 0:9
    
    DX.ComputeTBDDS; DX.ComputeTBDIS;
    DX.AddStatFluctTBDIS;
    
    TBDIS = DX.TBDIS;
    
    %TBDIS = D.TBDIS + (sqrt(D.TBDIS)).*randn(length(D.qU),1);
    qU = DX.qU + 0.5.*randn(length(DX.qU),1);
    
    TimeSec = DX.TimeSec + 10*rand - 5;
    qUfrac = DX.qUfrac;
    WGTS_CD_MolPerCm2 = DX.WGTS_CD_MolPerCm2 + 0.0001*DX.WGTS_CD_MolPerCm2*randn;
    WGTS_MolFrac_TT = DX.WGTS_MolFrac_TT + 0.001*DX.WGTS_MolFrac_TT*randn;
    WGTS_MolFrac_DT = DX.WGTS_MolFrac_DT + 0.001*DX.WGTS_MolFrac_DT*randn;
    WGTS_MolFrac_HT = DX.WGTS_MolFrac_HT + 0.001*DX.WGTS_MolFrac_HT*randn;
    
    save(['../../tritium-data/9999',num2str(ii),'.mat'],...
        'TBDIS','qU','TimeSec','qUfrac',...
        'WGTS_CD_MolPerCm2',...
        'WGTS_MolFrac_TT','WGTS_MolFrac_DT','WGTS_MolFrac_HT',...
        '-v7.3','-nocompression')
    
    TD = ['Run9999',num2str(ii)]; RunTime = TimeSec;
    
    save(['../../simulation/katrinsetup/TD_DataBank/Run9999',...
        num2str(ii),'.mat'],...
        'qU','qUfrac','RunTime','TD','-v7.3','-nocompression')
    
end