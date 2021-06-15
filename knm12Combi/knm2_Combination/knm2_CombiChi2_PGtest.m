
%% compatibility test: https://arxiv.org/abs/hep-ph/0304176v2
% calculate parameter goodness-of-fit
DataType = 'Real';
chi2 = 'chi2CMShape';
Knm2AnaFlag = 'Uniform';%MR-4';
nFit = 50;
%% load chi^2 profiles
d1 = LoadChi2Profile('DataSet','Knm1','DataType',DataType,'chi2',chi2,'AnaStr','Uniform','nFit',nFit,'mNuSqMin',-2.6,'mNuSqMax',1);
d2 = LoadChi2Profile('DataSet','Knm2','DataType',DataType,'chi2',chi2,'AnaStr',Knm2AnaFlag,'nFit',nFit,'mNuSqMin',-2.6,'mNuSqMax',1);

% load combi file (-> if not there: run knm2_CombiChi2.m)
ksumfile = [getenv('SamakPath'),sprintf('tritium-data/fit/Knm1/Chi2Profile/Uniform/Chi2ProfileCombi_%s_UniformScan_mNu_Knm1KNM2_UniformFPD_%s_FitParE0BkgNorm_nFit%.0f_min-2.6_max1.mat',DataType,chi2,nFit)];
dsum = importdata(ksumfile);
   
% asign variables
Chi21_min   = d1.ScanResults.BestFit.chi2;
Chi22_min   = d2.ScanResults.BestFit.chi2;
Chi2sum_min = dsum.BestFit.chi2;

mNuSq1_bf  = d1.ScanResults.BestFit.par;
mNuSq2_bf  = d2.ScanResults.BestFit.par;
mNuSq1_err = d1.ScanResults.BestFit.errMean;
mNuSq2_err = d2.ScanResults.BestFit.errMean;

%%

if strcmp(Knm2AnaFlag,'Uniform')
    N_x1 = 27;                  % subruns KNM1
    N_x2 = 28;                  % subruns KNM2
    dof1 = N_x1-4;
    dof2 = N_x2-4;
    N_tot = N_x1+N_x2;
    dof_tot = N_tot-7;
    dofHat = 4+4-7;            % not sure
elseif strcmp(Knm2AnaFlag,'MR-4')
    N_x1 = 4*27;                 % subruns KNM1
    N_x2 = 4*28;                 % subruns KNM2
    dof1 = N_x1-13;
    dof2 = N_x2-13;
    N_tot = N_x1+N_x2;
    dof_tot = N_tot-12;
    dofHat = 13+13-12;          % not sure
end
Chi2Hat_min = Chi2sum_min-Chi21_min-Chi22_min;
pHat = 1-chi2cdf(Chi2Hat_min,dofHat);
preg = 1-chi2cdf(Chi2sum_min,dof_tot);
p1 = 1-chi2cdf(Chi21_min,dof1);
p2 = 1-chi2cdf(Chi22_min,dof2);
%% result
fprintf('=================================================================== \n')
fprintf('parameter goodness-of-fit      p=%.2g (chi2minHat=%.1f, dofHat = %.2g)\n',pHat,Chi2Hat_min,dofHat);
fprintf('regular goodness-of-fit        p=%.2g (chi2min=%.1f  , dof    = %.2g)\n',preg,Chi2sum_min,dof_tot);
fprintf('regular goodness-of-fit KNM-1  p=%.2g (chi2min=%.1f  , dof    = %.2g)\n',p1,Chi21_min,dof1);
fprintf('regular goodness-of-fit KNM-2  p=%.2g (chi2min=%.1f  , dof    = %.2g)\n',p2,Chi22_min,dof2);

fprintf('=================================================================== \n')

%% look at chi2 
mNuSqCommon = 0.10;
KNM1SysBudget = 24;
KNM1Doppler   = 'OFF';

% KNM1
filecom = [getenv('SamakPath'),...
    sprintf('knm2ana/knm2_Combination/results/knm2_Combination_Chi2CommonPar_Knm1SysBudget%.0f_Knm1DE%s_mNuSqCommon%.3geV2.mat',...
    KNM1SysBudget,KNM1Doppler,mNuSqCommon)];
dcom = importdata(filecom);