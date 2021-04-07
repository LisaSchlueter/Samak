% dummy script to get relative uncertainty on variance
% uncertainties
nSamples = 5e4;
FSDShapeGS_RelErr = 0.04;
FSDShapeES_RelErr = 0.18;
FSD_NormRelErr = 0.01;

%% load fsd file
fsddir  = [getenv('SamakPath'),'inputs/FSD/'];
fsdfile = sprintf('%sFSD_KNM2_T2.txt',fsddir);
%fsdfile = sprintf('%sFSD_Saenz_T2mod.dat',fsddir);
d = importdata(fsdfile);
E    = d(:,1);
Prob = d(:,2);
EthIdx = find(E>=8,1);

%% randomization step I: normalization
NormTT_GS_i = sum(Prob(E<8));
NormTT_ES_i = sum(Prob(E>8));
NormTT_GS_Bias_rand = NormTT_GS_i.*randn(1,nSamples).*FSD_NormRelErr;
NormTT_GS_rand =  NormTT_GS_i+NormTT_GS_Bias_rand;
NormTT_ES_rand =  NormTT_ES_i+NormTT_GS_Bias_rand;

%% randomization step II: uncorrelated bin-to-bin fluctuations
Prob_rand_GS  = Prob(E<8).*(1+randn(sum(E<8),nSamples).*FSDShapeGS_RelErr);
Prob_rand_GS  = Prob_rand_GS.*NormTT_GS_rand./sum(Prob_rand_GS,1);
Prob_rand_ES  = Prob(E>8).*(1+randn(sum(E>8),nSamples).*FSDShapeES_RelErr);
Prob_rand_ES  = Prob_rand_ES.*NormTT_ES_rand./sum(Prob_rand_ES,1);
Prob_rand     = [Prob_rand_GS;Prob_rand_ES];

%% calculate variances etc. 
% boundedline(E,mean(Prob_rand,2),std(Prob_rand,0,2))
VarTotal_samples = zeros(5e3,1);
VarGS_samples    = zeros(5e3,1);
Var40_samples    = zeros(5e3,1);

for i=1:nSamples
VarGS_samples(i) = var(E(E<8),Prob_rand_GS(:,i));
Var40_samples(i) = var(E(E<=40),Prob_rand((E<=40),i));
VarTotal_samples(i) = var(E,Prob_rand(:,i));
end

% result
fprintf('Total var = %.1f eV^2 +-%.1f eV^2 (%.3g%%) \n',mean(VarTotal_samples),std(VarTotal_samples),1e2*std(VarTotal_samples)/VarTotal_mean)
fprintf('GS var = %.3f eV^2 +-%.2g eV^2  (%.2g%%) \n',mean(VarGS_samples),std(VarGS_samples),1e2*std(VarGS_samples)/VarGS_mean)
fprintf('40 eV var = %.1f eV^2 +-%.1f eV^2 (%.3f%%) \n',mean(Var40_samples),std(Var40_samples),1e2*std(Var40_samples)/Var40_mean)

