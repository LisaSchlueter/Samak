% 300 eV analysis
% calculate period-wise broadenings

nSamples = 1e4;%

savedir = [getenv('SamakPath'),'knm2ana/knm2_RateMonitoring/results/'];
savename = sprintf('%sknm2_RateMonitoring_GetCounts_%s.mat',savedir,'StackPixel');
load(savename);

%% period 1
R1 = RatesCorr(Idx_Rw1).*mean(TimeSec_ScanStep);
R1(abs(R1-mean(R1))>5*std(R1)) = [];
R1_rand =  R1(randi(numel(R1),numel(R1),nSamples));
R1_SigmaSqNorm = std(R1_rand).^2; 
R1_SigmaSqPoiss = mean(R1_rand);
R1_SigmaB_cps =  sqrt(R1_SigmaSqNorm-R1_SigmaSqPoiss)./mean(TimeSec_ScanStep);

%% period 2
R2 = RatesCorr(Idx_Rw2).*mean(TimeSec_ScanStep);
R2(abs(R2-mean(R2))>5*std(R2)) = [];
R2_rand =  R2(randi(numel(R2),numel(R2),nSamples));
R2_SigmaSqNorm = std(R2_rand).^2; 
R2_SigmaSqPoiss = mean(R2_rand);
R2_SigmaB_cps =  sqrt(R2_SigmaSqNorm-R2_SigmaSqPoiss)./mean(TimeSec_ScanStep);

%% period 3
R3 = RatesCorr(Idx_Rw3).*mean(TimeSec_ScanStep);
R3(abs(R3-mean(R3))>5*std(R3)) = [];
R3_rand =  R3(randi(numel(R3),numel(R3),nSamples));
R3_SigmaSqNorm = std(R3_rand).^2; 
R3_SigmaSqPoiss = mean(R3_rand);
R3_SigmaB_cps =  sqrt(R3_SigmaSqNorm-R3_SigmaSqPoiss)./mean(TimeSec_ScanStep);

%% meV equivalent
MeanRateAll = mean(RatesCorr);
qUSlope_absAll  = dqUCorr.par(1).*(mean(MeanRateAll)./mean(dqUCorr.Rate)); % cps/eV

R1_SigmaErr_mV = abs(1e3.*std(R1_SigmaB_cps)./qUSlope_absAll);
R2_SigmaErr_mV = abs(1e3.*std(R2_SigmaB_cps)./qUSlope_absAll);
R3_SigmaErr_mV = abs(1e3.*std(R3_SigmaB_cps)./qUSlope_absAll);

R1_Sigma_mV = abs(1e3.*mean(R1_SigmaB_cps)./qUSlope_absAll);
R2_Sigma_mV = abs(1e3.*mean(R2_SigmaB_cps)./qUSlope_absAll);
R3_Sigma_mV = abs(1e3.*mean(R3_SigmaB_cps)./qUSlope_absAll);

%% display
fprintf('Period 1: sigma = %.1f +- %.1f cps --> %.1f +- %.1f meV\n',...
    mean(R1_SigmaB_cps),std(R1_SigmaB_cps),R1_Sigma_mV,R1_SigmaErr_mV)
fprintf('Period 2: sigma = %.1f +- %.1f cps --> %.1f +- %.1f meV\n',...
    mean(R2_SigmaB_cps),std(R2_SigmaB_cps),R2_Sigma_mV,R2_SigmaErr_mV)
fprintf('Period 3: sigma = %.1f +- %.1f cps --> %.1f +- %.1f meV\n',...
    mean(R3_SigmaB_cps),std(R3_SigmaB_cps),R3_Sigma_mV,R3_SigmaErr_mV)



