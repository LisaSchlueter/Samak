

filedir = [getenv('SamakPath'),'tritium-data/fit/Knm2/Chi2Profile/Uniform/'];
filenameCommon = [filedir,'Chi2Profile_Real_Parallel_mNu_Knm2_UniformFPD_chi2CMShape_NP1.112_FitParmNuE0BkgNorm_nFit20_Min-1_Max1_KNM2_0p1eV_'];

ranges = [20:2:40,45:5:90];

GetFigure;
for i=1:numel(ranges)
dtmp = importdata(sprintf('%s%.0feVrange.mat',filenameCommon,ranges(i)));
% plot(dtmp.ParScan,dtmp.chi2min)
% hold on;

fprintf('Fit range %.0feV: mNuSq = %.3f +- %.3f eV^2 , E0 = %.2f +- %.2f eV\n',...
    ranges(i),dtmp.BestFit.par,dtmp.BestFit.errMean,dtmp.BestFit.E0+18573.7,dtmp.BestFit.E0Err);
end