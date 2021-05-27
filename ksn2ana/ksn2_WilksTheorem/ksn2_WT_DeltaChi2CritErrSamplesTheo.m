% estimate uncertainty on delta chi^2 crit. with bootstrapping
% using theoretical chi^2 function (2dof)
Mode = 'Theo';% use theoretical chi^2 function'Empirical';

nSamples = 2e3;
% re-sample from theoretical chi^2 distribution
x = linspace(0,30,2e3); % possible delta-chi^2 values
pdf =  chi2pdf(x,2);
cdf =  chi2cdf(x,2);
u = linspace(0,1,1e3);

nSamplesMax = [10:20:3000,3000];
nSteps = numel(nSamplesMax);

savefileErrCombi = sprintf('%sksn2ana/ksn2_WilksTheorem/results/Chi2CritConv/DeltaChi2CritErrConv_Theo_Max%.0fsamples_%.0fSteps.mat',...
    getenv('SamakPath'),max(nSamplesMax),nSteps);
if exist(savefileErrCombi,'file')
    load(savefileErrCombi);
else
    MeanChi2Crit_all = zeros(nSteps,1);
    StdChi2Crit_all   = zeros(nSteps,1);
    
    for i=1:nSteps
        progressbar(i./nSteps)
        savefileErr = sprintf('%sksn2ana/ksn2_WilksTheorem/results/Chi2CritConv/DeltaChi2CritErrConv_Theo_%.0fsamples.mat',...
            getenv('SamakPath'),nSamplesMax(i));
        
        if   exist(savefileErr,'file')
            d = importdata(savefileErr);
            MeanChi2Crit_all(i) =  d.MeanChi2Crit;
            StdChi2Crit_all(i) = d.StdChi2Crit;
        elseif ~exist(savefileErr,'file')
            
            
            DeltaChi2CDF_resample = zeros(nSamples,nSamplesMax(i));
            DeltaChi2Cr_reSample  = zeros(nSamples,1);
            DeltaChi2Cr_reSample2  = zeros(nSamples,1);
            
            chi2_delta_resample = interp1(cdf,x,rand(nSamplesMax(i),nSamples),'spline');
            
            for j=1:nSamples
                chi2_delta_tmp = chi2_delta_resample(:,j);
                DeltaChi2CDF_resample(j,:) = arrayfun(@(x) sum(chi2_delta_tmp<=x)./numel(chi2_delta_tmp),chi2_delta_tmp);
                DeltaChi2CDF_tmp = DeltaChi2CDF_resample(j,:);
                [DeltaChi2CDFquantile_tmp,ia] = unique(DeltaChi2CDF_tmp);
                DeltaChi2Cr_reSample(j) = interp1(DeltaChi2CDFquantile_tmp,chi2_delta_tmp(ia),0.95,'spline');
                DeltaChi2Cr_reSample2(j) = quantile(chi2_delta_tmp,0.95);
            end
            
            MeanChi2Crit = mean(DeltaChi2Cr_reSample);
            MedChi2Crit   = median(DeltaChi2Cr_reSample);
            StdChi2Crit   =std(DeltaChi2Cr_reSample);
            
            MeanChi2Crit_all(i) =  MeanChi2Crit;
            StdChi2Crit_all(i) = StdChi2Crit;
            fprintf('delta chi^2 crit = %.2f +- %.2f \n',mean(DeltaChi2Cr_reSample),std(DeltaChi2Cr_reSample))
            save(savefileErr,'MeanChi2Crit','StdChi2Crit');
        end
    end
    
    save(savefileErrCombi,'MeanChi2Crit_all','StdChi2Crit_all','nSamplesMax')
end
%%
%%
GetFigure;
plot(nSamplesMax,smooth(StdChi2Crit_all),'-','LineWidth',2)
xlabel('Samples')
ylabel(sprintf('\\Delta\\chi^2_{crit.}'))
PrettyFigureFormat;
savepltErr = sprintf('%sksn2ana/ksn2_WilksTheorem/plots/DeltaChi2CritErrConv_Theo_Max%.0fsamples.png',...
    getenv('SamakPath'),max(nSamplesMax));
leg = legend(sprintf('Uncertainty based on \\chi^2-function (2 dof)'));
PrettyLegendFormat(leg); leg.FontSize = get(gca,'FontSize');
print(savepltErr,'-dpng','-r350');

%%
Chi2CritTheo = GetDeltaChi2(0.95,2);

GetFigure;
plot(nSamplesMax,(MeanChi2Crit_all),'-','LineWidth',2)
hold on;
plot(nSamplesMax,Chi2CritTheo.*ones(nSteps,1),':k','LineWidth',2);
xlabel('Samples')
ylabel(sprintf('\\Delta\\chi^2'))
PrettyFigureFormat;

leg = legend(sprintf('Uncertainty based on \\chi^2-function (2 dof)'));
PrettyLegendFormat(leg); leg.FontSize = get(gca,'FontSize');
saveplt = sprintf('%sksn2ana/ksn2_WilksTheorem/plots/DeltaChi2CritConv_Theo_Max%.0fsamples.png',...
    getenv('SamakPath'),max(nSamplesMax));
print(saveplt,'-dpng','-r350');

