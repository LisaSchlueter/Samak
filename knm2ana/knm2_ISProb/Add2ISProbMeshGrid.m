function Add2ISProbMeshGrid(varargin)
p = inputParser;
p.addParameter('MACE_Bmax_T',4.23,@(x)isfloat(x));
p.addParameter('WGTS_B_T',2.52,@(x)isfloat(x));
p.addParameter('WGTS_CD_MolPerCm2',5*1e17,@(x)isfloat(x));

p.parse(varargin{:});
WGTS_CD_MolPerCm2 = p.Results.WGTS_CD_MolPerCm2;
MACE_Bmax_T       = p.Results.MACE_Bmax_T;
WGTS_B_T          = p.Results.WGTS_B_T;

savedir = [getenv('SamakPath'),'inputs/WGTSMACE/WGTS_ISProb/'];
savename = [savedir,'InitISProbMeshGrid.mat'];
dIniGrid = importdata(savename);

A = ref_FakeRun_KNM2_CD84_2hours('WGTS_CD_MolPerCm2',WGTS_CD_MolPerCm2,...
    'MACE_Bmax_T',MACE_Bmax_T,'WGTS_B_T',WGTS_B_T );

%%
maxEis = 1000;
IsProbBinStep = 10;
Eiscs = 18575+(-maxEis:IsProbBinStep:maxEis);


file_pis =  sprintf('%sIS_%.5g-molPercm2_Edep-Xsection-max%.0feV_Xstep%.1feV_%.0f-NIS_%.3g-Bmax_%.3g-Bs.mat',...
    savedir,A.WGTS_CD_MolPerCm2,maxEis,IsProbBinStep,A.NIS+1,A.MACE_Bmax_T,A.WGTS_B_T);
if ~exist(file_pis,'file')
    A.ComputeISProb('Energy',reshape(Eiscs,[1,1,numel(Eiscs)]));
end
d = importdata(file_pis);


[RhoDSigma,ia,ic] = unique([dIniGrid.RhoDSigma,A.WGTS_CD_MolPerCm2.*A.ISXsection(Eiscs)]);

ThetaFun = @(bmax,bs) asin(sqrt(bs./bmax));
Theta    = [dIniGrid.Theta,ThetaFun(A.MACE_Bmax_T,A.WGTS_B_T)];


ISProb = zeros(A.NIS+1,numel(Eiscs),size(dIniGrid.ISProb,3)+1);
ISProb(:,:,1:end-1)  = dIniGrid.ISProb;
ISProb(:,:,end)  = d.Pis_m;


%save(savename,'RhoDSigma','Theta','ISProb');
end