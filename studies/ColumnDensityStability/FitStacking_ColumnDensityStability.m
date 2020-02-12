

CDsigma_Rel = 0.0017;
nRuns = 500;

% load MC
savedir = [getenv('SamakPath'),'studies/ColumnDensityStability/results/'];
MakeDir(savedir);
savename = [savedir,sprintf('ColumnDensityStability_RelStd%.2fprcnt_nRuns%.0f.mat',CDsigma_Rel*100,nRuns)];

d = importdata(savename);
   
% calculate model but for longer measurement time
InitFile = @ref_FakeRun_FinalKATRIN_CD100_2hours;

M = RunAnalysis('DataType','Fake',...
                 'RunNr',1,...
                 'FakeInitFile',InitFile,...
                 'TwinBias_Time',2*60*60*nRuns,...
                 'fixPar','mNu E0 Norm Bkg',...
                 'PixList',1:124);
M.Fit;
mNuSq_ref = M.FitResult.par(1);

%%
M.RunData.TBDIS = sum(d.TBDIS,2);
M.Fit;
mNuSq = M.FitResult.par(1);