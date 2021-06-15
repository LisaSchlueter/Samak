

savedir = [getenv('SamakPath'),'tritium-data/FC/DeltaChi2LookupTable/'];
% savename1 =[savedir,'AsimovDeltaChi2_mNuSq0eV2_KNM1_chi2CMShape_40bE0_freeParmNuE0BkgNorm_300samples.mat'];
% %savename2 = [savedir,'AsimovDeltaChi2_mNuSq0eV2_KNM2_Prompt_chi2CMShape_SysBudget40_41bE0_freeParmNuE0BkgNorm_300samples.mat'];
% d = importdata(savename1);
% mNuSq = d.mNuSq;
% Chi2 = d.Chi2True;

savename3 = [getenv('SamakPath'),'tritium-data/fit/Knm1/Chi2Profile/Uniform/Chi2Profile_Real_UniformScan_mNu_Knm1_UniformFPD_chi2CMShape_SysBudget24_NP1.064_FitParE0BkgNorm_nFit200_min-5_max5.mat'];
d = importdata(savename3);
[mNuSq,Idx] = unique(reshape(d.ScanResults.ParScan,numel(d.ScanResults.ParScan),1));
Chi2  = reshape(d.ScanResults.chi2min,numel(d.ScanResults.ParScan),1);
Chi2 = Chi2(Idx);

Prob = exp(-0.5.*Chi2)./simpsons(mNuSq,exp(-0.5.*Chi2));
ProbCDF = GetCDF(mNuSq,Prob);

close all;
GetFigure
plot(mNuSq,Prob,'LineWidth',2)
hold on;
plot(mNuSq,ProbCDF);
grid on
xlim([-5 5])
PrettyFigureFormat;
xlabel(sprintf('{\\itm}^2 (eV^2)'))


Write2Txt('filename','KNM1_Data_Chi2Profile','Format','dat','variable',[mNuSq,Chi2,Prob],'variableName','mNuSq Chi2 Prob','nCol',3);