% neutrino mass bias as a function of column density stability 
% 1000 fake runs with different column density
% column density is randomly drawn from Gaussian distribution
% stack runs and get neutrino mass bias ("shift")
% settings: 2 hour runs, 100% column density
% compute Fake runs
%%
CDsigma_Rel = 0.0017;
nRuns = 500;

% label 
savedir = [getenv('SamakPath'),'studies/ColumnDensityStability/results/'];
MakeDir(savedir);
savename = [savedir,sprintf('ColumnDensityStability_RelStd%.2fprcnt_nRuns%.0f.mat',CDsigma_Rel*100,nRuns)];
   
%%
WGTS_CD_MolPerCm2 = 5*1e17.*(1+CDsigma_Rel.*randn(nRuns,1));

% calculate model
InitFile = @ref_FakeRun_FinalKATRIN_CD100_2hours;
D = InitFile();

%%
TBDIS = zeros(D.nqU,nRuns);

for i=1:numel(WGTS_CD_MolPerCm2)
    progressbar(i/numel(WGTS_CD_MolPerCm2))
    D.WGTS_CD_MolPerCm2 = WGTS_CD_MolPerCm2(i);
    D.AdjustRF;
    D.ComputeTBDDS;
    D.ComputeTBDIS;
    TBDIS(:,i) = D.TBDIS;
end

save(savename,'TBDIS','WGTS_CD_MolPerCm2','InitFile');
