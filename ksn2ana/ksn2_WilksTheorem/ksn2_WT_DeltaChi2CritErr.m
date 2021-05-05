% estimate uncertainty on delta chi^2 crit. with bootstrapping
Hypothesis = 'H0';
InterpMode = 'lin';
SavePlt = 'ON';
MergeNew = 'ON';
RmDuplicates = 'ON';

switch Hypothesis
    case 'H0'
        randMC =[1001:1260,1294:1300,1349:1500];%11:1e3;
        Twin_sin2T4 = 0;
        Twin_mNu4Sq = 0;
        chi2 = 'chi2CMShape';
        randMC_new  = [1:96,110:200,406:500];
    case 'H1' 
        randMC = 1:1e3;
        Twin_sin2T4 = 0.0240;
        Twin_mNu4Sq = 92.7;
        chi2 = 'chi2Stat';
        MergeNew = 'OFF'; % nothing new
end

if strcmp(MergeNew,'ON')
    MergeStr = sprintf('_MergeNew%.0f',numel(randMC_new));
      NrandMC = numel(randMC)+numel(randMC_new);
else
    MergeStr = '';
    NrandMC = numel(randMC);
end

savedir = [getenv('SamakPath'),'ksn2ana/ksn2_WilksTheorem/results/'];

if Twin_sin2T4==0 && Twin_mNu4Sq==0
    savefile = sprintf('%sksn2_WilksTheorem_NullHypothesis_Interp%s_%.0fsamples%s_RmDouble%s.mat',...
        savedir,InterpMode,numel(randMC),MergeStr,RmDuplicates);
else
    savefile = sprintf('%sksn2_WilksTheorem_mNu4Sq-%.1feV2_sin2T4-%.3g_Interp%s_%.0fsamples%s_RmDouble%s.mat',...
        savedir,Twin_mNu4Sq,Twin_sin2T4,InterpMode,numel(randMC),MergeStr,RmDuplicates);
end

if exist(savefile,'file')
    load(savefile);
    fprintf('load file from %s \n',savefile);
else
     fprintf('file does not exist: %s \n',savefile);
     return
end

%%
% mean: calculate 95 quantile: interpolation
PlotDeltaChi2 = sort(chi2_delta);
DeltaChi2CDF = arrayfun(@(x) sum(PlotDeltaChi2<=x)./numel(PlotDeltaChi2),PlotDeltaChi2);
[DeltaChi2CDFquantile,ia] = unique(DeltaChi2CDF);
DeltaChi2Cr = interp1(DeltaChi2CDFquantile,PlotDeltaChi2(ia),0.95,'lin');%quantile(PlotDeltaChi2,0.95);% PlotDeltaChi2(find(abs(DeltaChi2CDF-0.95)==min(abs(DeltaChi2CDF-0.95)),1));

%% bootstrapping error: re-sampling with double counting
nSamples = 5e3;
tmp = repmat(chi2_delta,1,nSamples);
chi2_delta_resample = tmp(randi(numel(chi2_delta),numel(chi2_delta),nSamples));
DeltaChi2Cr_reSample = zeros(nSamples,1);
for i=1:nSamples
    chi2_delta_tmp = sort(chi2_delta_resample(:,i));
    DeltaChi2CDF_tmp = arrayfun(@(x) sum(chi2_delta_tmp<=x)./numel(chi2_delta_tmp),chi2_delta_tmp);
    [DeltaChi2CDFquantile_tmp,ia] = unique(DeltaChi2CDF_tmp);
    DeltaChi2Cr_reSample(i) = interp1(DeltaChi2CDFquantile_tmp,chi2_delta_tmp(ia),0.95,'lin');
end

fprintf('delta chi^2 crit = %.2f +- %.2f \n',mean(DeltaChi2Cr_reSample),std(DeltaChi2Cr_reSample))
