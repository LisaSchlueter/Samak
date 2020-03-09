function [RhoDSigma,Theta,ISProb] = Compute_InitISProbMeshGrid
savedir = [getenv('SamakPath'),'inputs/WGTSMACE/WGTS_ISProb/'];
savename = [savedir,'InitISProbMeshGrid.mat'];


if exist(savename,'file')
    load(savename,'RhoDSigma','Theta','ISProb')
else
    A = ref_FakeRun_KNM2_CD84_2hours;
    
    %%
    maxEis = 1000;
    IsProbBinStep = 2;
    Eiscs = 18575+(-maxEis:IsProbBinStep:maxEis);
    Bmax_Min = A.MACE_Bmax_T.*0.90;
    Bmax_Max = A.MACE_Bmax_T.*1.10;
    BmaxSamples = linspace(Bmax_Min,Bmax_Max,20);
    
    ISProb = zeros(A.NIS+1,numel(Eiscs),numel(BmaxSamples));
    
    for i=1:numel(BmaxSamples)
        progressbar(i/numel(BmaxSamples));
        file_pis =  sprintf('%sIS_%.5g-molPercm2_Edep-Xsection-max%.0feV_Xstep%.1feV_%.0f-NIS_%.3g-Bmax_%.3g-Bs.mat',...
            savedir,A.WGTS_CD_MolPerCm2,maxEis,IsProbBinStep,A.NIS+1,BmaxSamples(i),A.WGTS_B_T);
        if ~exist(file_pis,'file')
            A.MACE_Bmax_T = BmaxSamples(i);
            A.ComputeISProb('Energy',reshape(Eiscs,[1,1,numel(Eiscs)]));
        end
        
        d = importdata(file_pis);
        ISProb(:,:,i) = d.Pis_m(:,:);
    end
    
    
    RhoDSigma = A.WGTS_CD_MolPerCm2.*A.ISXsection(Eiscs);
    ThetaFun = @(bmax,bs)  asin(sqrt(bs./bmax));
    Theta    = ThetaFun(BmaxSamples,A.WGTS_B_T);
    save(savename,'RhoDSigma','Theta','ISProb');
end
end