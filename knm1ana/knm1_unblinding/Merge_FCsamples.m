
mNuSq_all = [-0.98];%0.1,0.15,0.2:0.1:0.7,0.9,1.2,1.35,1.5];
nSamples = 5000;
for i=1:numel(mNuSq_all)
dir = [getenv('SamakPath'),'tritium-data/fit/TwinMCKnm1/'];

str = sprintf('TwinMC%.3geV2mNuSq_KNM1_chi2CMShape_40bE0_fixPar567891011_%.0fsamples.mat',mNuSq_all(i),nSamples);
fprintf('Loading %s from file \n',str);
d1 = importdata([dir,str]);
d2 = importdata([dir,sprintf('TwinMC%.3geV2mNuSq_KNM1_chi2CMShape_40bE0_fixPar567891011_%.0fsamples.mat',mNuSq_all(i),nSamples+1)]);
%d2 = importdata([dir,'TwinMC0eV2mNuSq_KNM1_chi2Stat_40bE0_fixPar567891011_1501samples.mat']);

par = [d1.par,d2.par(:,1:nSamples)];
err = [d1.err,d2.err(:,1:nSamples)];
chi2min = [d1.chi2min;d2.chi2min(1:nSamples)];
dof = d1.dof;
mNuSq = d1.mNuSq;
%Chi2True = [d1.Chi2True';d2.Chi2True(1:1500)']';
%Chi2Best = [d1.Chi2Best';d2.Chi2Best(1:1500)']';
%DeltaChi2 = [d1.DeltaChi2';d2.DeltaChi2(1:1500)'];


savename = [dir,sprintf('TwinMC%.2geV2mNuSq_KNM1_chi2CMShape_40bE0_fixPar567891011_%.0fsamples.mat',mNuSq_all(i),nSamples*2)];
save(savename,'par','err','chi2min','dof','mNuSq');%,'Chi2True','Chi2Best','DeltaChi2');
disp(savename)
end
%%


