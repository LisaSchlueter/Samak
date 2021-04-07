%% 1. likelihood from fit on randomized twins
savedir = [getenv('SamakPath'),'knm1ana/knm1_Likelihood/results/'];
savename = sprintf('%sknm1_MarginalLikelihood.mat',savedir);
mNuSq_fit = importdata(savename);
%h1 = histogram(mNuSq);

%% likelihood from chi2 profile
savename3 = [getenv('SamakPath'),'tritium-data/fit/Knm1/Chi2Profile/Uniform/Chi2Profile_Real_UniformScan_mNu_Knm1_UniformFPD_chi2CMShape_SysBudget24_NP1.064_FitParE0BkgNorm_nFit200_SibilleFull_min-5_max5.mat'];
d = importdata(savename3);
[mNuSq,Idx] = unique(reshape(d.ScanResults.ParScan,numel(d.ScanResults.ParScan),1));
Chi2  = reshape(d.ScanResults.chi2min,numel(d.ScanResults.ParScan),1);
Chi2 = Chi2(Idx);
Prob = exp(-0.5.*Chi2)./simpsons(mNuSq,exp(-0.5.*Chi2));
ProbCDF = GetCDF(mNuSq,Prob);

%% plot
GetFigure;
h1 = histogram(mNuSq_fit,'Normalization','pdf','FaceColor',rgb('SkyBlue'),'FaceAlphA',1,'EdgeColor',rgb('PowderBlue'));
hold on;
p1 = plot(mNuSq,Prob,'LineWidth',2.5);
PrettyFigureFormat;
xlabel(sprintf('{\\itm}_\\nu^2 (eV^2)'));
xlim([-5,3]);
ylabel('Probability density (eV^{-2})');
leg = legend([h1,p1],sprintf('1000 MC fits ({\\itm}^2_{true}= -0.97 eV^2)'),sprintf('Profile Likelihood from \\chi^2 profile (data)'));
PrettyLegendFormat(leg);

