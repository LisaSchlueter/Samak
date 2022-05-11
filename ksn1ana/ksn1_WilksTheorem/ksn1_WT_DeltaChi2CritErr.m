range = 40;%
SavePlt = 'OFF';
chi2Str = 'chi2CMShape';
InterpMode = 'lin';
savedir = [getenv('SamakPath'),'ksn1ana/ksn1_WilksTheorem/results/'];
MakeDir(savedir);
savename = sprintf('%sksn1_WilksTheorem_%.0frange_%s_%s.mat',savedir,range,chi2Str,InterpMode);

if exist(savename,'file')
    d = importdata(savename);
    fprintf('load %s\n',savename);
else
    fprintf(2,'file not found, run ksn1_WT_MergeFiles.m \n');
    return
end

chi2_delta = d.chi2min_null-d.chi2min_bf;
%% DeltaChi2 distribution (best fit - null chi2)
dof = 2;
chi2_delta(chi2_delta<0) = 0;
PlotDeltaChi2 = sort(chi2_delta);%unique(chi2_delta);%
DeltaChi2CDF = arrayfun(@(x) sum(PlotDeltaChi2<=x)./numel(PlotDeltaChi2),PlotDeltaChi2);

% calculate 95 quantile: interpolation
[DeltaChi2CDFquantile,ia] = unique(DeltaChi2CDF);
DeltaChi2CrApprox = interp1(DeltaChi2CDFquantile,PlotDeltaChi2(ia),0.95,'lin');%quantile(PlotDeltaChi2,0.95);% PlotDeltaChi2(find(abs(DeltaChi2CDF-0.95)==min(abs(DeltaChi2CDF-0.95)),1));
x = linspace(0,max(PlotDeltaChi2),1e3);
DeltaChi2CDFTheo = chi2cdf(x,dof);
xInter = linspace(0,10,1e2);
DeltaChi2CrTheo = interp1(chi2cdf(xInter,dof),xInter,0.95,'spline');

%% mean: calculate 95 quantile: interpolation
chi2_delta(chi2_delta<0) = 0;
PlotDeltaChi2 = sort(chi2_delta);
DeltaChi2CDF = arrayfun(@(x) sum(PlotDeltaChi2<=x)./numel(PlotDeltaChi2),PlotDeltaChi2);
[DeltaChi2CDFquantile,ia] = unique(DeltaChi2CDF);
DeltaChi2Cr = interp1(DeltaChi2CDFquantile,PlotDeltaChi2(ia),0.95,'spline');%quantile(PlotDeltaChi2,0.95);% PlotDeltaChi2(find(abs(DeltaChi2CDF-0.95)==min(abs(DeltaChi2CDF-0.95)),1));

if strcmp(Mode,'Empirical')
    %% bootstrapping error: re-sampling with double counting
    nSamples = 5e3;
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
% GetFigure;
% histogram(DeltaChi2Cr_reSample,'FaceAlpha',0.2)
% hold on
% histogram(DeltaChi2Cr_reSample2,'FaceAlpha',0.2)
% hold off;
% legend('Interp','Quantile')