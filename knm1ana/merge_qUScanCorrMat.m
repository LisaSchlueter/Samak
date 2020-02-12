RunList = 'KNM1';
firstPoint = 2;
lastPoint = 20;
chi2 = 'chi2CMShape';
nSamples1 = 300;
nSamples2 = 201;

savedir = [getenv('SamakPath'),'knm1ana/knm1_qUScan/results/'];
savefile1 = [savedir,sprintf('CorrMat_qUScan_%s_%.0f-%.0f_%s_%.0f.mat',RunList,firstPoint,lastPoint,chi2,nSamples1)];
savefile2 = [savedir,sprintf('CorrMat_qUScan_%s_%.0f-%.0f_%s_%.0f.mat',RunList,firstPoint,lastPoint,chi2,nSamples2)];

d1 = importdata(savefile1);
d2 = importdata(savefile2);

chi2qU = [d1.chi2qU,d2.chi2qU];
dofqU  = [d1.dofqU,d2.dofqU];

parqU = zeros(size(d1.parqU,1),size(d1.parqU,2),nSamples2+nSamples1);
parqU(:,:,1:nSamples1) = d1.parqU;
parqU(:,:,nSamples1+1:end) = d2.parqU;

errqU = zeros(size(d1.errqU,1),size(d1.errqU,2),nSamples2+nSamples1);
errqU(:,:,1:nSamples1) = d1.errqU;
errqU(:,:,nSamples1+1:end) = d2.errqU;

CorrMatmNuSq = corr(squeeze(parqU(1,:,:))');
CorrMatE0 = corr(squeeze(parqU(2,:,:))');
CorrMatB = corr(squeeze(parqU(3,:,:))');
CorrMatN = corr(squeeze(parqU(4,:,:))');

savefile3 = [savedir,sprintf('CorrMat_qUScan_%s_%.0f-%.0f_%s_%.0f.mat',RunList,firstPoint,lastPoint,chi2,nSamples2+nSamples1)];
save(savefile3,'parqU','errqU','chi2qU','dofqU','CorrMatmNuSq','CorrMatE0','CorrMatB','CorrMatN');

