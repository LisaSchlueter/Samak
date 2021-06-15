% marginal likelihood
% fit randomized MC data
% script that does the fits: knm1_FitRandData.m
%% settings
Plot = 'OFF';
Save = 'ON';
TwinmNuSq = -0.97;
savedirRandFit = [getenv('SamakPath'),'knm1ana/knm1_FitRandData/results/'];
savedir = [getenv('SamakPath'),'knm1ana/knm1_Likelihood/results/'];
savefileCommon = sprintf('%sknm1FitRandData_Twin_chi2CMShape_SysBudget24_NP1.064_mNuE0NormBkg_40eV_Sibille0p5eV_TwinmNuSq%.3geV2',savedirRandFit,TwinmNuSq);
d1 = importdata([savefileCommon,'_1000fit.mat']);
d2 = importdata([savefileCommon,'_1001fit.mat']);
mNuSq = [d1.FitPar(1,:),d2.FitPar(1,:)];
MakeDir(savedir);

fprintf('Mean   = %.3f \n',mean(mNuSq));
fprintf('Median = %.3f \n',median(mNuSq));

if strcmp(Plot,'ON')
    h1 = histogram(mNuSq);
    PrettyFigureFormat;
    xlabel('m^2');
end

if strcmp(Save,'ON')
    savename = sprintf('%sknm1_MarginalLikelihood.mat',savedir);
    save(savename,'mNuSq');
end
