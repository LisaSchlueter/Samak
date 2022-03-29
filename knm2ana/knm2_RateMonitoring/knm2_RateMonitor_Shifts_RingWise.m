% 300 eV analysis
% calculate shifts from corrected data

% load data
savedir = [getenv('SamakPath'),'knm2ana/knm2_RateMonitoring/results/'];
savename = sprintf('%sknm2_RateMonitoring_GetCounts_%s%s.mat',savedir,'Ring','Full');
load(savename);

% mean rates per period
MeanRate1 = mean(RatesCorr(:,Idx_Rw1),2);
MeanRate2 = mean(RatesCorr(:,Idx_Rw2),2);
MeanRate3 = mean(RatesCorr(:,Idx_Rw3),2);
MeanRateAll = mean(RatesCorr,2);

% error of the mean
EoMeanRate1 = std(RatesCorr(:,Idx_Rw1),0,2)./sqrt(sum(Idx_Rw1));
EoMeanRate2 = std(RatesCorr(:,Idx_Rw2),0,2)./sqrt(sum(Idx_Rw2));
EoMeanRate3 = std(RatesCorr(:,Idx_Rw3),0,2)./sqrt(sum(Idx_Rw3));
EoMeanRateAll = std(RatesCorr,0,2)./sqrt(numel(Idx_Rw1));

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
fprintf('                               Period1                    Period 2:               Period3 \n');
for i=1:numel(MeVDiff1)
    
    fprintf('Ring %.2g (%.0f) meV diff:   %.1f +- %.1f meV             %.1f +- %.1f meV          %.1f +- %.1f meV \n',...
        i,numel(MeVDiff1),MeVDiff1(i),MeVDiffErr1(i),MeVDiff2(i),MeVDiffErr2(i),MeVDiff3(i),MeVDiffErr3(i));
end


fprintf('Max radial diff:           %.1f meV                       %.1f meV                  %.1f meV \n',...
    max(MeVDiff1)-min(MeVDiff1), max(MeVDiff2)-min(MeVDiff2), max(MeVDiff3)-min(MeVDiff3));




