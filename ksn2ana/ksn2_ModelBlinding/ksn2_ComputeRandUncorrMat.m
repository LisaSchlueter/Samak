% compute and save additional covariance matrix for ksn-2 model blinding test
% uncorrelated (diagonal only)
% relative uncertainties -> fractional
% mu = mean additional uncertainty
% sigma = std of addition uncertainty (random over subruns)
mu    = 0.0005;
sigma = 0.001;
nSubrun = 38;
AddCovMatFrac = eye(nSubrun).*(sigma.*randn(nSubrun,1)+mu).^2;

savedir  = [getenv('SamakPath'),'ksn2ana/ksn2_ModelBlinding/results/'];
savefile = sprintf('%sksn2_RandUncorrFracMat_mu%.3g_std%.3g.mat',savedir,mu,sigma);
MakeDir(savedir);
save(savefile,'AddCovMatFrac');