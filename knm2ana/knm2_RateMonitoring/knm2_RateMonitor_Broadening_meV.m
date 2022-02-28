% conver broadening, calculated in knm2_RateMonitor_Broadening.m in meV
% equivalent
savename = sprintf('%sknm2_RateMonitoring_GetCounts_%s.mat',savedir,'StackPixel');
load(savename);
savedir = [getenv('SamakPath'),'knm2ana/knm2_RateMonitoring/results/'];

SigmaSqNorm  = zeros(3,1);
SigmaSqPoiss = zeros(3,1);

for i=1:3
    tmpname = sprintf('%sknm2_RateMonitor_Broadening_RW%.0f.mat',savedir,i);
    d = importdata(tmpname);
    SigmaSqNorm(i) = d.SigmaSqNorm_t;
    SigmaSqPoiss(i) = d.SigmaSqPoiss_t;
end

%% broadening in cps
SigmaB_cps =  sqrt(SigmaSqNorm-SigmaSqPoiss)./mean(TimeSec_ScanStep);

% conversion to mV
MeanRateAll = mean(RatesCorr);
qUSlope_absAll  = dqUCorr.par(1).*(mean(MeanRateAll)./mean(dqUCorr.Rate)); % cps/eV
SigmaB_mV = abs(1e3.*SigmaB_cps./qUSlope_absAll);




