% 300 eV analysis
% calculate shifts from corrected data

% load data
savedir = [getenv('SamakPath'),'knm2ana/knm2_RateMonitoring/results/'];
savename = sprintf('%sknm2_RateMonitoring_GetCounts_%s.mat',savedir,'StackPixel');
load(savename);

% mean rates per period
MeanRate1 = mean(RatesCorr(Idx_Rw1));
MeanRate2 = mean(RatesCorr(Idx_Rw2));
MeanRate3 = mean(RatesCorr(Idx_Rw3));
MeanRateAll = mean(RatesCorr);

% error of the mean
EoMeanRate1 = std(RatesCorr(Idx_Rw1))./sqrt(sum(Idx_Rw1));
EoMeanRate2 = std(RatesCorr(Idx_Rw2))./sqrt(sum(Idx_Rw2));
EoMeanRate3 = std(RatesCorr(Idx_Rw3))./sqrt(sum(Idx_Rw3));
EoMeanRateAll = std(RatesCorr)./sqrt(numel(Idx_Rw1));

% rate difference
RateDiff1 = (MeanRate1-MeanRateAll);
RateDiff2 = (MeanRate2-MeanRateAll);
RateDiff3 = (MeanRate3-MeanRateAll);

% uncertainty on rate difference (simple Gaussian error prop. of error-of-mean)
RateDiffErr1 = EoMeanRate1;%sqrt(EoMeanRate1^2+EoMeanRateAll^2);
RateDiffErr2 = EoMeanRate2;%sqrt(EoMeanRate2^2+EoMeanRateAll^2);
RateDiffErr3 = EoMeanRate3;%sqrt(EoMeanRate3^2+EoMeanRateAll^2);

% conversion to meV
qUSlope_absAll  = dqUCorr.par(1).*(mean(MeanRateAll)./mean(dqUCorr.Rate)); % cps/eV

MeVDiff1 = 1e3.*RateDiff1./qUSlope_absAll;
MeVDiff2 = 1e3.*RateDiff2./qUSlope_absAll;
MeVDiff3 = 1e3.*RateDiff3./qUSlope_absAll;

MeVDiffErr1 = abs(1e3.*RateDiffErr1./qUSlope_absAll);
MeVDiffErr2 = abs(1e3.*RateDiffErr2./qUSlope_absAll);
MeVDiffErr3 = abs(1e3.*RateDiffErr3./qUSlope_absAll);
%% rates:
fprintf('Period1: Rates diff =   %.1f +- %.1f cps -> meV diff:  %.1f +- %.1f meV\n',...
    RateDiff1,RateDiffErr1,MeVDiff1,MeVDiffErr1);
fprintf('Period2: Rates diff =  %.1f +- %.1f cps -> meV diff: %.1f +- %.1f meV\n',...
    RateDiff2,RateDiffErr2,MeVDiff2,MeVDiffErr2);
fprintf('Period3: Rates diff = %.1f +- %.1f cps -> meV diff:  %.1f +- %.1f meV\n',...
    RateDiff3,RateDiffErr3,MeVDiff3,MeVDiffErr3);





