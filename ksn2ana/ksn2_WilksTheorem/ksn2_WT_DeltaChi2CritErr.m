% estimate uncertainty on delta chi^2 crit. with bootstrapping
Hypothesis = 'H0';
SavePlt = 'ON';
MergeNew = 'ON';
RmDuplicates = 'ON';

Mode = 'Empirical';%'Theo';% use theoretical chi^2 function'Empirical';

switch Hypothesis
    case 'H0'
        randMC =[1001:1260,1294:1300,1349:1500];%11:1e3;
        Twin_sin2T4 = 0;
        Twin_mNu4Sq = 0;
        chi2 = 'chi2CMShape';
        randMC_new  = 1:1250;
        InterpMode = 'lin';
    case 'H1' 
       randMC = 1:1500;
        Twin_sin2T4 = 0.0240;
        Twin_mNu4Sq = 92.7;
        chi2 = 'chi2CMShape';
        MergeNew = 'OFF'; % nothing new
        InterpMode = 'Mix';
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
    savefile = sprintf('%sksn2_WilksTheorem_mNu4Sq-%.1feV2_sin2T4-%.3g_Interp%s_%.0fsamples.mat',...
        savedir,Twin_mNu4Sq,Twin_sin2T4,InterpMode,numel(randMC));
end

if exist(savefile,'file')
    load(savefile);
    fprintf('load file from %s \n',savefile);
else
     fprintf('file does not exist: %s \n',savefile);
     return
end

%% mean: calculate 95 quantile: interpolation
chi2_delta(chi2_delta<0) = 0;
PlotDeltaChi2 = sort(chi2_delta);
DeltaChi2CDF = arrayfun(@(x) sum(PlotDeltaChi2<=x)./numel(PlotDeltaChi2),PlotDeltaChi2);
[DeltaChi2CDFquantile,ia] = unique(DeltaChi2CDF);
DeltaChi2Cr = interp1(DeltaChi2CDFquantile,PlotDeltaChi2(ia),0.95,'spline');%quantile(PlotDeltaChi2,0.95);% PlotDeltaChi2(find(abs(DeltaChi2CDF-0.95)==min(abs(DeltaChi2CDF-0.95)),1));

if strcmp(Mode,'Empirical')
    %% bootstrapping error: re-sampling with double counting
    nSamples = 10e3;
    tmp                   = repmat(chi2_delta,1,nSamples);
    chi2_delta_resample   = tmp(randi(numel(chi2_delta),numel(chi2_delta),nSamples))+rand(numel(chi2_delta),nSamples).*1e-5; %
    chi2_delta_resample   = sort(chi2_delta_resample,1);
    DeltaChi2Cr_reSample  = zeros(nSamples,1);
    DeltaChi2Cr_reSample2  = zeros(nSamples,1);
    DeltaChi2CDF_resample = zeros(nSamples,numel(chi2_delta));
    
    for j=1:nSamples
        %  chi2_delta_resample(:,j) = sort(chi2_delta_resample(:,j));
        chi2_delta_tmp = chi2_delta_resample(:,j);
        DeltaChi2CDF_resample(j,:) = arrayfun(@(x) sum(chi2_delta_tmp<=x)./numel(chi2_delta_tmp),chi2_delta_tmp);
        DeltaChi2CDF_tmp = DeltaChi2CDF_resample(j,:);
        [DeltaChi2CDFquantile_tmp,ia] = unique(DeltaChi2CDF_tmp);
        DeltaChi2Cr_reSample(j) = interp1(DeltaChi2CDFquantile_tmp,chi2_delta_tmp(ia),0.95,'spline');
        DeltaChi2Cr_reSample2(j) = quantile(chi2_delta_tmp,0.95);
    end
    MeanChi2Crit = mean(DeltaChi2Cr_reSample);
    MedChi2Crit = median(DeltaChi2Cr_reSample);
    StdChi2Crit =std(DeltaChi2Cr_reSample);
else
    nSamples = 2e3;
    % re-sample from theoretical chi^2 distribution
    x = linspace(0,30,2e3); % possible delta-chi^2 values
    pdf =  chi2pdf(x,2);
    cdf =  chi2cdf(x,2);
    u = linspace(0,1,1e3);
    chi2_delta_resample = interp1(cdf,x,rand(numel(chi2_delta),nSamples),'spline');
    DeltaChi2CDF_resample = zeros(nSamples,numel(chi2_delta));
     
    DeltaChi2Cr_reSample  = zeros(nSamples,1);
    DeltaChi2Cr_reSample2  = zeros(nSamples,1);
    for j=1:nSamples
        chi2_delta_tmp = chi2_delta_resample(:,j);
        DeltaChi2CDF_resample(j,:) = arrayfun(@(x) sum(chi2_delta_tmp<=x)./numel(chi2_delta_tmp),chi2_delta_tmp);
        DeltaChi2CDF_tmp = DeltaChi2CDF_resample(j,:);
        [DeltaChi2CDFquantile_tmp,ia] = unique(DeltaChi2CDF_tmp);
        DeltaChi2Cr_reSample(j) = interp1(DeltaChi2CDFquantile_tmp,chi2_delta_tmp(ia),0.95,'spline');
        DeltaChi2Cr_reSample2(j) = quantile(chi2_delta_tmp,0.95);
    end
    
    MeanChi2Crit = mean(DeltaChi2Cr_reSample);
    MedChi2Crit = median(DeltaChi2Cr_reSample);
    StdChi2Crit =std(DeltaChi2Cr_reSample);  
%     GetFigure
%     h1 = histogram(deltachi2_samples(:,1),'Normalization','probability');
%     hold on;
%     plot(x,pdf.*h1.BinWidth)
end
%%
fprintf('delta chi^2 crit = %.2f (mean) %.2f (median) +- %.2f \n',MeanChi2Crit,MedChi2Crit,StdChi2Crit)
GetFigure;
histogram(DeltaChi2Cr_reSample,'FaceAlpha',0.2)
hold on
histogram(DeltaChi2Cr_reSample2,'FaceAlpha',0.2)
hold off;
legend('Interp','Quantile')




