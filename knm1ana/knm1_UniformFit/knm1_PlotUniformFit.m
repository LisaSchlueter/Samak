% plot knm1 final uniform fit

%% load result: stat & syst
% knm1: prl config stat & syst.

savedir = [getenv('SamakPath'),'knm12Combi/knm1ReAnalysis/results/'];
f_prl_Tot = sprintf('%sknm1_UniformFit_ReAna_Real_chi2CMShape_SysBudget24_NP1.064_mNuE0NormBkg_40eV_Sibille0p5eV_KatrinT2_AngTFOFF.mat',savedir);
d_prl_tot = importdata(f_prl_Tot);
mNuSqTot(1) = d_prl_tot.FitResult.par(1);
mNuSqErrTot(1) = 0.5*(d_prl_tot.FitResult.errPos(1)-d_prl_tot.FitResult.errNeg(1));

%%
d_prl_tot.Real.PlotFit('FitResultsFlag','OFF','YLimRes',[-2.2 2.7],'saveplot','pdf');