% find out impact of gas composition uncertainty
kappa_Mean = 4;
epsT_Mean = 0.975;
kappa_Std = 0.05; % kappa = WGTS_MolFrac_HT./WGTS_MolFrac_DT;
epsT_Std = 0.001; % WGTS_espT = WGTS_MolFrac_TT+0.5.*(WGTS_MolFrac_HT+WGTS_MolFrac_DT);
nSamples = 1e3;


savedir = [getenv('SamakPath'),'knmxana/knmx_GasComposition/results/'];
savefile = [savedir,sprintf('knmx_GasComposition_epsT_%.3f_err%.3g_kappa_%.3g_err_%.3g_%.0fSamples.mat',...
    epsT_Mean,epsT_Std,kappa_Mean,kappa_Std,nSamples)];
if exist(savefile,'file')
load(savefile)
else
%%
%% sample kappa and tritium purity
kappa_s = kappa_Mean.*(1+kappa_Std.*randn(nSamples,1));
epsT_s = epsT_Mean.*(1+epsT_Std.*randn(nSamples,1));

% translate into molecular fractions
WGTS_MolFrac_TT = 2.*epsT_s-1;
WGTS_MolFrac_HT = (1-WGTS_MolFrac_TT)./(1+1./kappa_s);
WGTS_MolFrac_DT = (1-WGTS_MolFrac_TT-WGTS_MolFrac_HT);

%% plot
close all
GetFigure
scatter(WGTS_MolFrac_TT,WGTS_MolFrac_HT)
hold on;
scatter(WGTS_MolFrac_TT,WGTS_MolFrac_DT)

A = ref_KNMX_Final_2h_GasCombo;
A.ComputeTBDDS;A.ComputeTBDIS;
TBDIS_i = A.TBDIS;
%
TBDIS_s = zeros(39,nSamples);
A_s = ref_KNMX_Final_2h_GasCombo;
for i=1:nSamples
    A_s.WGTS_MolFrac_TT = WGTS_MolFrac_TT(i);
    A_s.WGTS_MolFrac_HT = WGTS_MolFrac_HT(i);
    A_s.WGTS_MolFrac_DT = WGTS_MolFrac_DT(i);
    A_s.LoadFSD;
    A_s.ComputeTBDDS;A_s.ComputeTBDIS;
    TBDIS_s(:,i) = A_s.TBDIS;
end

%%
R = RunAnalysis('RunNr',1,...
    'DataType','Fake','FakeInitFile',...
    @ref_KNMX_Final_2h_1000runs_GasCombo,'chi2','chi2Stat',...
    'ELossFlag','KatrinT2A20','FSDFlag','KNM2_0p5eV',...
    'FixPar','mNu E0 Norm Bkg',...
    'exclDataStart',13,'DopplerEffectFlag','FSD');
R.RunData.TBDIS = 1000.*TBDIS_i;
R.Fit;
FitResults_i = R.FitResult;

R.RunData.TBDIS = sum(TBDIS_s,2);
R.Fit;
FitResults_s = R.FitResult;

save('WGTS_MolFrac_TT','WGTS_MolFrac_HT','WGTS_MolFrac_DT','kappa_s','epsT_s',...
    'TBDIS_i','TBDIS_s','FitResults_i','FitResults_s');
end

