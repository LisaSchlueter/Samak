% KNM2 LARA Data
savedir  = [getenv('SamakPath'),'knm1ana/knm1_GasComposition/results/'];
savefile = sprintf('%sknm1_LARAData.mat',savedir);

if exist(savefile,'file') 
    d = importdata(savefile);
else
    M = MultiRunAnalysis('RunList','KNM1','NonPoissonScaleFactor',1);
    WGTS_MolFrac_TT = M.SingleRunData.WGTS_MolFrac_TT_SubRun;
    WGTS_MolFrac_HT = M.SingleRunData.WGTS_MolFrac_HT_SubRun;
    WGTS_MolFrac_DT = M.SingleRunData.WGTS_MolFrac_DT_SubRun;
    
    nBin = M.ModelObj.nqU.*M.nRuns;
    WGTS_MolFrac = [reshape(WGTS_MolFrac_TT,nBin,1),...
        reshape(WGTS_MolFrac_HT,nBin,1),...
        reshape(WGTS_MolFrac_DT,nBin,1)];
    CorrMat = corrcoef(WGTS_MolFrac);
    MakeDir(savedir);
    save(savefile,'CorrMat','WGTS_MolFrac_TT','WGTS_MolFrac_HT','WGTS_MolFrac_DT','WGTS_MolFrac');
end





