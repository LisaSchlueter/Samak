% find out impact of gas composition uncertainty
% fit model: mc truth
% simulation: true-+trueness shift
kappa_Mean = 4;
epsT_Mean = 0.975;
kappa_Std = 0.033; % kappa = WGTS_MolFrac_HT./WGTS_MolFrac_DT;
epsT_Std = 0.0016; % WGTS_espT = WGTS_MolFrac_TT+0.5.*(WGTS_MolFrac_HT+WGTS_MolFrac_DT);
nSamples = 1e3;


savedir = [getenv('SamakPath'),'knmxana/knmx_GasComposition/results/'];
savefile = [savedir,sprintf('knmx_GasComposition_epsT_%.3f_err%.3g_kappa_%.3g_err_%.3g_Shift.mat',...
    epsT_Mean,epsT_Std,kappa_Mean,kappa_Std)];
if exist(savefile,'file')
load(savefile)
else
%%
%% sample kappa and tritium purity
kappa_s = kappa_Mean.*(1+[kappa_Std,-kappa_Std]);
epsT_s = epsT_Mean.*(1+[epsT_Std,-epsT_Std]);

% translate into molecular fractions
WGTS_MolFrac_TT = 2.*epsT_s-1;
WGTS_MolFrac_HT = (1-WGTS_MolFrac_TT)./(1+1./kappa_s);
WGTS_MolFrac_DT = (1-WGTS_MolFrac_TT-WGTS_MolFrac_HT);


A = ref_KNMX_Final_2h_1000runs_GasCombo;
A.ComputeTBDDS;A.ComputeTBDIS;
TBDIS_i = A.TBDIS;
%
TBDIS_sim = zeros(39,2);
for i=1:2
   A_sim = ref_KNMX_Final_2h_1000runs_GasCombo('WGTS_MolFrac_TT',WGTS_MolFrac_TT(i),...
    'WGTS_MolFrac_HT',WGTS_MolFrac_HT(i),...
    'WGTS_MolFrac_DT',WGTS_MolFrac_DT(i));
    A_sim.ComputeTBDDS;A_sim.ComputeTBDIS;
    TBDIS_sim(:,i) = A_sim.TBDIS;
end

%%
R = RunAnalysis('RunNr',1,...
    'DataType','Fake','FakeInitFile',...
    @ref_KNMX_Final_2h_1000runs_GasCombo,'chi2','chi2Stat',...
    'ELossFlag','KatrinT2A20','FSDFlag','KNM2_0p5eV',...
    'FixPar','mNu E0 Norm Bkg',...
    'exclDataStart',13,'DopplerEffectFlag','FSD');
R.RunData.TBDIS = TBDIS_i;
R.Fit;
FitResults_i = R.FitResult;
%%
R.RunData.TBDIS = TBDIS_sim(:,1);
R.Fit;
FitResults_sim1 = R.FitResult;

R.RunData.TBDIS = TBDIS_sim(:,2);
R.Fit;
FitResults_sim2 = R.FitResult;

save('WGTS_MolFrac_TT','WGTS_MolFrac_HT','WGTS_MolFrac_DT','kappa_s','epsT_s',...
    'TBDIS_i','TBDIS_sim','FitResults_i','FitResults_sim1','FitResults_sim2');
end

