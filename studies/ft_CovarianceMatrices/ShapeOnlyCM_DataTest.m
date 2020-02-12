% Script to Test Shape Only Analysis with Simulated Data
RunList = 'StackCD100all';
exclDataStart = 7;
nSamples = 1000;

R = MultiRunAnalysis('RunList',RunList,'fixPar','1 5 6','exclDataStart',exclDataStart,'chi2','chi2CM');
R.Fit;
FitResultsCM = R.FitResult;
R.chi2 = 'chi2CMShape';
R.Fit;
FitResultsCMShape = R.FitResult;

%% plot
figure(12);
subplot(2,2,1);
imagesc(R.FitCM)
colorbar;

subplot(2,2,2);
imagesc(R.FitCMShape);
colorbar;

