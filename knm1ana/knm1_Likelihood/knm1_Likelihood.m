

%% 1. likelihood from fit on randomized twins
nFit = 1000;
RecomputeFlag = 'OFF';
Plots        = 'OFF';
DataType     = 'Twin';
range        = 40;
chi2         = 'chi2CMShape';
freePar      = 'mNu E0 Norm Bkg';
FSDFlag      = 'Sibille0p5eV';
NonPoissonScaleFactor = 1.064;
SysBudget    = 24;
savedir = [getenv('SamakPath'),'knm1ana/knm1_FitRandData/results/'];
TwinmNuSq = -0.97;
savefile = sprintf('%sknm1FitRandData_%s_%s_NP%.4g_%s_%.0feV_%s_TwinmNuSq%.3geV2_%.0ffit.mat',...
    savedir,DataType,chi2,NonPoissonScaleFactor,strrep(freePar,' ',''),range,FSDFlag,TwinmNuSq,nFit);
if strcmp(chi2,'chi2CMShape')
    savefile = strrep(savefile,chi2,sprintf('%s_SysBudget%.0f',chi2,SysBudget));
end
dfits = importdata(savefile); 
mNuSq_fit = dfits.FitPar(1,:);
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

