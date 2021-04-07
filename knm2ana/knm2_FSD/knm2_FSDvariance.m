% look at FSD and calculate variance
fsddir  = [getenv('SamakPath'),'inputs/FSD/'];
fsdfile = sprintf('%sFSD_KNM2_T2.txt',fsddir);
d = importdata(fsdfile);
E = d(:,1);
Prob = d(:,2);
fprintf('----------------T2-----------------\n')
fprintf('Ground State variance: %.3f eV^2 \n', var(E(E<=8),Prob(E<=8))); % variance ground state
fprintf('40 eV variance:        %.3f eV^2 \n',var(E(E<=40),Prob(E<=40))); % variance ground state
fprintf('Total variance:        %.3f eV^2 \n', var(E,Prob)); % variance ground state


fsddir  = [getenv('SamakPath'),'inputs/FSD/'];
fsdfile = sprintf('%sFSD_KNM2_HT.txt',fsddir);
d = importdata(fsdfile);
E = d(:,1);
Prob = d(:,2);
fprintf('----------------HT-----------------\n')
fprintf('Ground State variance: %.3f eV^2 \n', var(E(E<=8),Prob(E<=8))); % variance ground state
fprintf('40 eV variance:        %.3f eV^2 \n',var(E(E<=40),Prob(E<=40))); % variance ground state
fprintf('Total variance:        %.3f eV^2 \n', var(E,Prob)); % variance ground state

fsddir  = [getenv('SamakPath'),'inputs/FSD/'];
fsdfile = sprintf('%sFSD_KNM2_DT.txt',fsddir);
d = importdata(fsdfile);
E = d(:,1);
Prob = d(:,2);
fprintf('----------------DT-----------------\n')
fprintf('Ground State variance: %.3f eV^2 \n', var(E(E<=8),Prob(E<=8))); % variance ground state
fprintf('40 eV variance:        %.3f eV^2 \n',var(E(E<=40),Prob(E<=40))); % variance ground state
fprintf('Total variance:        %.3f eV^2 \n', var(E,Prob)); % variance ground state


%%
d = importdata('FSD_KNM2_DT-HT-TT-CovMat_5000Trials_KNM2Prompt_0.01NormErr_0.04GS_0.18ES_ShapeErr.mat');
E   =  d.obj.StudyObject.TTexE;
Prob = squeeze(d.TT_P_norm);
Prob_mean = mean(Prob,2);

VarTotal_mean = var(E,Prob_mean);
VarGS_mean = var(E(E<8),Prob_mean(E<8));
Var40_mean = var(E(E<=40),Prob_mean(E<=40));

VarTotal_samples = zeros(5e3,1);
VarGS_samples    = zeros(5e3,1);
Var40_samples    = zeros(5e3,1);

for i=1:5e3
VarTotal_samples(i) = var(E,Prob(:,i));
VarGS_samples(i) = var(E(E<8),Prob((E<8),i));
Var40_samples(i) = var(E(E<=40),Prob((E<=40),i));
end
%%
fprintf('Total var = %.1f eV^2 +-%.1f eV^2 (%.2g%%) \n',VarTotal_mean,std(VarTotal_samples),1e2*std(VarTotal_samples)/VarTotal_mean)
fprintf('GS var = %.1f eV^2 +-%.2g eV^2(%.2g%%) \n',VarGS_mean,std(VarGS_samples),1e2*std(VarGS_samples)/VarGS_mean)
fprintf('40 eV var = %.1f eV^2 +-%.1f eV^2 (%.2g%%) \n',Var40_mean,std(Var40_samples),1e2*std(Var40_samples)/Var40_mean)


