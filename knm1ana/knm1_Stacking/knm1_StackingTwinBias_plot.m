

% fit to different kind of twins average response function for KNM1
RunList = 'KNM1';
[Twin, ~, ~, ~, ~, ~,~] = ComputeLoadTwinObjects('RunList',RunList);
savedir = [getenv('SamakPath'),'knm1ana/knm1Twins/results/'];

FitRange40 = round(-Twin.ModelObj.qU(14)+18573.7); %for labeling
save_file = [savedir,'FitResult_Twins_all',RunList,sprintf('%.0feVrange',FitRange40),'.mat'];
d40 = importdata(save_file);
mNuSq40     = d40.mNuSq;
mNuSqErr40  = d40.mNuSqErr;
E040        = d40.E0;
E0Err40     = d40.E0Err;
chi2min40   = d40.chi2min;
dof40 = d40.dof;

FitRange30 = round(-Twin.ModelObj.qU(17)+18573.7); %for labeling
save_file = [savedir,'FitResult_Twins_all',RunList,sprintf('%.0feVrange',FitRange30),'.mat'];
d30 = importdata(save_file);
mNuSq30     = d30.mNuSq;
mNuSqErr30  = d30.mNuSqErr;
E030        = d30.E0;
E0Err30     = d30.E0Err;
chi2min30   = d30.chi2min;
dof30 = d30.dof;