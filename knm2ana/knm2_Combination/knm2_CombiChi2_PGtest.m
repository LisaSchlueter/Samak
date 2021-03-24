
%% compatibility test: https://arxiv.org/abs/hep-ph/0304176v2
% calculate parameter goodness-of-fit
DataType = 'Real';
chi2 = 'chi2CMShape';
Knm2AnaFlag = 'Uniform';%MR-4';

%% load chi^2 profiles
k1file = [getenv('SamakPath'),sprintf('tritium-data/fit/Knm1/Chi2Profile/Uniform/Chi2Profile_%s_UniformScan_mNu_Knm1_UniformFPD_%s_NP1.064_FitParE0BkgNorm_nFit20_min-2_max1.mat',DataType,chi2)];
d1 = importdata(k1file); fprintf('load knm1: %s \n',k1file);

if strcmp(Knm2AnaFlag,'Uniform')
k2file = [getenv('SamakPath'),sprintf('tritium-data/fit/Knm2/Chi2Profile/Uniform/Chi2Profile_%s_UniformScan_mNu_Knm2_UniformFPD_%s_NP1.112_FitParE0BkgNorm_nFit20_min-2_max1.mat',DataType,chi2)];
elseif strcmp(Knm2AnaFlag,'MR-4')
k2file = [getenv('SamakPath'),sprintf('tritium-data/fit/Knm2/Chi2Profile/Ring_Full/Chi2Profile_%s_UniformScan_mNu_Knm2_Ring_FullFPD_%s_NP1.112_FitParE0BkgNormqU_nFit20_min-2_max1.mat',DataType,chi2)];   
end
d2 = importdata(k2file); fprintf('load knm2: %s \n',k2file)

ksumfile = [getenv('SamakPath'),sprintf('tritium-data/fit/Knm1/Chi2Profile/Uniform/Chi2ProfileCombi_%s_UniformScan_mNu_Knm1KNM2_UniformFPD_%s_FitParE0BkgNorm_nFit20_min-2_max1.mat',DataType,chi2)];
dsum = importdata(ksumfile);
   
% asign variables
Chi21_min = d1.BestFit.chi2;
Chi22_min = d2.BestFit.chi2;
Chi2sum_min = dsum.BestFit.chi2;

mNuSq1_bf  = d1.BestFit.mNuSq;
mNuSq2_bf  = d2.BestFit.mNuSq;
mNuSq1_err = d1.BestFit.mNuSqErr;
mNuSq2_err = d2.BestFit.mNuSqErr;

%%

if strcmp(Knm2AnaFlag,'Uniform')
    N_x1 = 27;                  % subruns KNM1
    N_x2 = 28;                  % subruns KNM2
    N_tot = N_x1+N_x2;
    dof_tot = N_tot-7;
    dofHat = 4+4-7;            % not sure
elseif strcmp(Knm2AnaFlag,'MR-4')
    N_x1 = 4*27;                 % subruns KNM1
    N_x2 = 4*28;                 % subruns KNM2
    N_tot = N_x1+N_x2;
    dof_tot = N_tot-12;
    dofHat = 13+13-12;          % not sure
end
Chi2Hat_min = Chi2sum_min-Chi21_min-Chi22_min;
pHat = 1-chi2cdf(Chi2Hat_min,dofHat);
preg = 1-chi2cdf(Chi2sum_min,dof_tot);
%% result
fprintf('=================================================================== \n')
fprintf('parameter goodness-of-fit p=%.2g (chi2minHat=%.1f, dofHat = %.2g)\n',pHat,Chi2Hat_min,dofHat);
fprintf('regular goodness-of-fit   p=%.2g (chi2min=%.1f  , dof    = %.2g)\n',preg,Chi2sum_min,dof_tot);
fprintf('=================================================================== \n')

%% look at chi2 
filecom = [getenv('SamakPath'),sprintf('knm2ana/knm2_Combination/results/knm2_Combination_Chi2CommonPar.mat')];
dcom = importdata(filecom);