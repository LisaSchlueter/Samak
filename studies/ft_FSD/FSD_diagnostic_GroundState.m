% Study to test FSD systematics treatment after Hamish Robertson's suggestion
% Introduce bin-to-bin uncorrelated fluctuation such that the variance of
% the FSD will have a variation (std) of 1% on ground state
Isotopologue = 'TT';
PlotFontSize = 18;

% Configuration
FSDNorm_RelErr = 0.01;
FSDShapeGS_RelErr = [0.2, 0.3];%0.01, 0.02, 0.025 0.03, 0.035, 0.04, 0.05, 0.06];%[0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04];
FSDShapeES_RelErr = zeros(numel(FSDShapeGS_RelErr),1);

nTrials = 10000;
state = 'ground state';
%save_name = 'blalala';
save_name = sprintf('./results/%s_FSDdiagnostic_GroundState_%.0fSamples.mat',Isotopologue,nTrials);
if exist(save_name,'file')
    load(save_name);
    close;
else
% Set up Model 
A = MultiRunAnalysis('RunList','StackCD100all');

switch Isotopologue
    case 'DT'
        nProb   = numel(A.ModelObj.DTexP);
        Energy  = A.ModelObj.DTexE;
        FSDCalc = A.ModelObj.DTFSD;
        
    case 'TT'
        nProb   = numel(A.ModelObj.TTexP);
        Energy  = A.ModelObj.TTexE;
        FSDCalc = A.ModelObj.TTFSD;
    case 'HT'
        nProb   = numel(A.ModelObj.HTexP);
        Energy  = A.ModelObj.HTexE;
        FSDCalc = A.ModelObj.HTFSD;
end
GSESLimit =A.ModelObj.TTGSTh; %for all same bin

% Init Variables
FSD_P_norm = zeros(nProb,nTrials,numel(FSDShapeGS_RelErr)); % Final State Probabilities after norm (happends inside ComputeTBDDS)
FSD_P      = zeros(nProb,nTrials,numel(FSDShapeGS_RelErr)); % Final State Probabilities before norm
NormGS = zeros(nTrials,numel(FSDShapeGS_RelErr));                            % Normalization Factor Ground State
NormES = zeros(nTrials,numel(FSDShapeGS_RelErr));                            % Normalization Factor Excited State
NormBias = zeros(nTrials,numel(FSDShapeGS_RelErr));

%Get theortical FSD
A.ModelObj.ComputeTBDDS;
%Get theortical FSD
A.ModelObj.ComputeTBDDS;
switch Isotopologue
    case 'DT'
        FSD_P_theory = [A.ModelObj.DTexP_G*A.ModelObj.DTNormGS_i,A.ModelObj.DTexP_E*A.ModelObj.DTNormES_i]';
        FSD_Norm_i = A.ModelObj.DTNormGS_i;
    case 'TT'
        FSD_P_theory = [A.ModelObj.TTexP_G*A.ModelObj.TTNormGS_i,A.ModelObj.TTexP_E*A.ModelObj.TTNormES_i]';
        FSD_Norm_i = A.ModelObj.TTNormGS_i;
    case 'HT'
        FSD_P_theory = [A.ModelObj.HTexP_G*A.ModelObj.HTNormGS_i,A.ModelObj.HTexP_E*A.ModelObj.HTNormES_i]';
        FSD_Norm_i = A.ModelObj.HTNormGS_i;
end

for ii=1:numel(FSDShapeGS_RelErr)

% Do Variation
NormBias(:,ii) = randn(nTrials,1).*FSD_Norm_i*FSDNorm_RelErr;

switch Isotopologue
    case 'DT'
        FSD_P(1:A.ModelObj.DTGSTh,:,ii)     = A.ModelObj.DTexP(1:A.ModelObj.DTGSTh)'.*(1+randn(numel(A.ModelObj.DTexP(1:A.ModelObj.DTGSTh)),nTrials).*FSDShapeGS_RelErr(ii));
        FSD_P(A.ModelObj.DTGSTh+1:end,:,ii) = A.ModelObj.DTexP(A.ModelObj.DTGSTh+1:end)'.*(1+randn(numel(A.ModelObj.DTexP(A.ModelObj.DTGSTh+1:end)),nTrials).*FSDShapeES_RelErr(ii));
    case 'TT'
        FSD_P(1:A.ModelObj.TTGSTh,:,ii)     = A.ModelObj.TTexP(1:A.ModelObj.TTGSTh)'.*(1+randn(numel(A.ModelObj.TTexP(1:A.ModelObj.TTGSTh)),nTrials).*FSDShapeGS_RelErr(ii));
        FSD_P(A.ModelObj.TTGSTh+1:end,:,ii) = A.ModelObj.TTexP(A.ModelObj.TTGSTh+1:end)'.*(1+randn(numel(A.ModelObj.TTexP(A.ModelObj.TTGSTh+1:end)),nTrials).*FSDShapeES_RelErr(ii));
    case 'HT'
         FSD_P(1:A.ModelObj.HTGSTh,:,ii)     = A.ModelObj.HTexP(1:A.ModelObj.HTGSTh)'.*(1+randn(numel(A.ModelObj.HTexP(1:A.ModelObj.HTGSTh)),nTrials).*FSDShapeGS_RelErr(ii));
         FSD_P(A.ModelObj.HTGSTh+1:end,:,ii) = A.ModelObj.HTexP(A.ModelObj.HTGSTh+1:end)'.*(1+randn(numel(A.ModelObj.HTexP(A.ModelObj.HTGSTh+1:end)),nTrials).*FSDShapeES_RelErr(ii));
end
FSD_P(FSD_P<0)=0; %set all negative probabilities to 0

% Compute Model with varied parameters
for i=1:nTrials
    switch Isotopologue
        case 'DT'
            A.ModelObj.DTexP = FSD_P(:,i,ii)';
            A.ModelObj.ComputeTBDDS('DTGS_bias',NormBias(i,ii),'DTES_bias', -NormBias(i,ii));
            FSD_P_norm(:,i,ii) = [A.ModelObj.DTexP_G*A.ModelObj.DTNormGS,A.ModelObj.DTexP_E*A.ModelObj.DTNormES;]';
        case 'TT'
            A.ModelObj.TTexP = FSD_P(:,i,ii)';
            A.ModelObj.ComputeTBDDS('TTGS_bias',NormBias(i,ii),'TTES_bias', -NormBias(i,ii));
            FSD_P_norm(:,i,ii) = [A.ModelObj.TTexP_G*A.ModelObj.TTNormGS,A.ModelObj.TTexP_E*A.ModelObj.TTNormES;]';
        case 'HT'
            A.ModelObj.HTexP = FSD_P(:,i,ii)';
            A.ModelObj.ComputeTBDDS('HTGS_bias',NormBias(i,ii),'HTES_bias', -NormBias(i,ii));
            FSD_P_norm(:,i,ii) = [A.ModelObj.HTexP_G*A.ModelObj.HTNormGS,A.ModelObj.HTexP_E*A.ModelObj.HTNormES;]';
    end
end
end

% variance distribution of samples widht different sampling uncertainties
FSDmean = zeros(nTrials,numel(FSDShapeGS_RelErr));  % weighted mean energy   (eV)
FSDvar  = zeros(nTrials,numel(FSDShapeGS_RelErr));  % weighted mean variance of FSD (eV^2)
FSDvar_GS  = zeros(nTrials,numel(FSDShapeGS_RelErr));  % weighted mean variance of FSD (eV^2) for ground state
meanFSDvar = zeros(numel(FSDShapeGS_RelErr),1);     % mean variance of FSD variance(=width)
meanFSDvar_GS = zeros(numel(FSDShapeGS_RelErr),1);  
varFSDvar  = zeros(numel(FSDShapeGS_RelErr),1);     % variance of FSD variance(=width)
varFSDvar_GS = zeros(numel(FSDShapeGS_RelErr),1);

for ii=1:numel(FSDShapeGS_RelErr)
for i=1:nTrials
    FSDmean(i,ii) = wmean(Energy,FSD_P_norm(:,i,ii)');
    FSDvar(i,ii) = var(Energy,FSD_P_norm(:,i,ii)');
    FSDvar_GS(i,ii) = var(Energy(1:A.ModelObj.DTGSTh), FSD_P_norm(1:A.ModelObj.DTGSTh,i,ii)');
end
meanFSDvar(ii)   = mean(FSDvar(:,ii));
meanFSDvar_GS(ii) = mean(FSDvar_GS(:,ii));
varFSDvar(ii)    = var(FSDvar(:,ii));
varFSDvar_GS(ii) = var(FSDvar_GS(:,ii));
end

FSD_P_norm_std  = squeeze(std(permute(FSD_P_norm,[2,1,3]))); % standard deviation of samles for plot
FSD_P_norm_mean = squeeze(mean(permute(FSD_P_norm,[2,1,3])));% mean of samles for plot
%save(save_name)
end
%% FSD distribution for different sampling uncertainties
mytitle = sprintf('%u DT Sample Distributions with \nnormalization uncertainty %u %% and different bin-to-bin uncorrelated uncertainties (ground state)',nTrials, FSDNorm_RelErr*100);
for ii=1:numel(FSDShapeGS_RelErr)
    IndexFSDRelErr = ii;
    fig55 = figure('Name',mytitle,'Renderer','opengl');
    set(fig55, 'Units', 'normalized', 'Position', [0.9, 0.9, 1, 0.82]);
    [l,p] = boundedline(Energy,FSD_P_norm_mean(:,IndexFSDRelErr)*100 ,FSD_P_norm_std(:,IndexFSDRelErr)*100);
    l.LineWidth = 1.5; l.Color = rgb('CadetBlue');
    p.FaceColor = rgb('CadetBlue'); p.FaceAlpha = 0.5;
    set(gca,'XScale','log')
    xlabel('excitation energy (eV)');
    ylabel('probability (%)');
    %title(sprintf('%u DT Sample Distributions with \nnormalization uncertainty %u %%',nTrials, FSDNorm_RelErr*100));
    title(sprintf('normalization uncertainty 1%%, bin-to-bin uncorrelated uncertainty (GS) %.1f%%', FSDShapeGS_RelErr(IndexFSDRelErr)*100));
    sample_leg = sprintf('%u Sample Distributions',nTrials);
    theo_leg = sprintf('%s Theory FSD (%s)',Isotopologue,FSDCalc);
    legend([l p],theo_leg, sample_leg,'Location','northeast');
    legend boxoff;
    PrettyFigureFormat;
    xlim([min(Energy) max(Energy)]);%max(A.ModelObj.DTGSTh-25)]);% max(Energy)]);
    ylim([min(FSD_P_theory*100) 8.2]);%max(FSD_P_norm_mean(:,IndexFSDRelErr)*100)+max(FSD_P_norm_std(:,IndexFSDRelErr))*100]);
    set(gca,'FontSize',PlotFontSize);
      
save_dir = sprintf('./plots/FSD_diagnostic/');
save_name = sprintf('SampleFSD_GroundStates_%.2gGSerr_%.0fSamples',FSDShapeGS_RelErr(IndexFSDRelErr),nTrials);
export_fig([save_dir,'png/',save_name,'.png']);
export_fig([save_dir,'fig/',save_name,'.fig']);
publish_figurePDF(fig55,[save_dir,'pdf/',save_name,'.pdf'])
close;
end
%% Variance distribution 
close
fig67 = figure('Name','Variances','Renderer','opengl');
IndexFSDRelErr = 5;
set(fig67, 'Units', 'normalized', 'Position', [0.9, 0.9, 0.75, 0.7]);
    h1 = nhist(FSDvar_GS(:,IndexFSDRelErr),'color',rgb('CadetBlue'),'linewidth',0.01);
    sample_title = sprintf('bin-to-bin uncorrelated uncertainty GS %.1f%% , normalization uncertainty %.0f %%', FSDShapeGS_RelErr(IndexFSDRelErr)*100,FSDNorm_RelErr(1)*100);
    sample_leg = sprintf('\\sigma_{FSD}=%.1f %%',sqrt(varFSDvar_GS(IndexFSDRelErr))./meanFSDvar_GS(IndexFSDRelErr)*100);
    l1 = legend(sample_leg,'Location','northwest');
   legend boxoff
    maintitle=  sprintf('\\mu=%.3f eV^2, \\sigma=%.3feV^2',meanFSDvar_GS(IndexFSDRelErr),sqrt(varFSDvar_GS(IndexFSDRelErr)));
    xlabel('variance \sigma^2 (eV²)'); ylabel('samples');
    title(sample_title);
    PrettyFigureFormat;
    set(gca,'FontSize',PlotFontSize);
    l1.FontSize = PlotFontSize+2;

a=annotation('textbox', [0 0.2 1 0], ...
    'String', maintitle, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
a.FontSize=PlotFontSize;a.FontWeight='bold';    
    
save_dir = sprintf('./plots/FSD_diagnostic/');
save_name = sprintf('Variance_GroundStates_%.1f_%.0fSamples',FSDShapeGS_RelErr(IndexFSDRelErr)*100,nTrials);
export_fig([save_dir,'png/',save_name,'.png']);
export_fig([save_dir,'fig/',save_name,'.fig']);
publish_figurePDF(fig67,[save_dir,'pdf/',save_name,'.pdf']);
%close;
%% Variance in dependence of GS Error
fig681 = figure(681);
set(fig681, 'Units', 'normalized', 'Position', [0.9, 0.9, 0.7, 0.7],'Renderer','opengl');
plot(FSDShapeGS_RelErr*100,sqrt(varFSDvar_GS)./meanFSDvar_GS*100,'--o','MarkerSize',8,'LineWidth',5,'Color',rgb('FireBrick'),'MarkerFaceColor',rgb('FireBrick'));
hold on;
plot(FSDShapeGS_RelErr*100,ones(numel(FSDShapeGS_RelErr),1),'--k','LineWidth',2);
xlabel(sprintf('\\sigma_{Fluct} (%%)'));
ylabel('\sigma_{FSD} (%)');
PrettyFigureFormat;
ylim([min(sqrt(varFSDvar_GS)./meanFSDvar_GS*100) max(sqrt(varFSDvar_GS)./meanFSDvar_GS*100)])
title(sprintf('%s  Ground State Sampling (%u Samples) \nNormalization Uncertainty %u %%',Isotopologue,nTrials, FSDNorm_RelErr*100));
set(gca,'FontSize',PlotFontSize);
grid on;

save_dir = sprintf('./plots/FSD_diagnostic/');
save_name = sprintf('%s_Variance_OverviewRatio_GroundStates_%.0fSamples',Isotopologue,nTrials);
export_fig([save_dir,'png/',save_name,'.png']);
export_fig([save_dir,'fig/',save_name,'.fig']);
publish_figurePDF(fig681,[save_dir,'pdf/',save_name,'.pdf']);
%close;
%% std of width of FSD depending in sampling uncertainty Overview
fig68 = figure('Renderer','opengl');
set(fig68, 'Units', 'normalized', 'Position', [0.9, 0.9, 0.5, 1]);
s1 = subplot(3,1,1);
plot(FSDShapeGS_RelErr*100,sqrt(varFSDvar_GS),'--o','MarkerSize',8,'LineWidth',3,'Color',rgb('FireBrick'),'MarkerFaceColor',rgb('FireBrick'));
xlabel(sprintf('FSD bin-to-bin uncorrelated uncertainty (%s)',state));
ylabel('\sigma of FSD variance (eV^2)')
PrettyFigureFormat;
set(gca,'FontSize',10);
ylim([min(sqrt(varFSDvar_GS)) max(sqrt(varFSDvar_GS))]);
title(sprintf('%s FSD Ground State Sampling (%u Samples) \nnormalization uncertainty %u %%',Isotopologue,nTrials, FSDNorm_RelErr*100));
grid on;

s2 = subplot(3,1,3);
plot(FSDShapeGS_RelErr*100,sqrt(varFSDvar_GS)./meanFSDvar_GS*100,'--o','MarkerSize',8,'LineWidth',3,'Color',rgb('FireBrick'),'MarkerFaceColor',rgb('FireBrick'));
hold on;
plot(FSDShapeGS_RelErr*100,ones(numel(FSDShapeGS_RelErr),1),'--k');
xlabel(sprintf('FSD bin-to-bin uncorrelated uncertainty (%s)',state));
ylabel('\sigma of FSD variance / mean variance (%)');
PrettyFigureFormat;
set(gca,'FontSize',10);
ylim([min(sqrt(varFSDvar_GS)./meanFSDvar_GS*100) max(sqrt(varFSDvar_GS)./meanFSDvar_GS*100)])
grid on;

s3 = subplot(3,1,2);
plot(FSDShapeGS_RelErr*100,meanFSDvar_GS*100,'--o','MarkerSize',8,'LineWidth',3,'Color',rgb('FireBrick'),'MarkerFaceColor',rgb('FireBrick'));
xlabel(sprintf('FSD bin-to-bin uncorrelated uncertainty (%s)',state));
ylabel('mean FSD variance (eV^2)')
PrettyFigureFormat;
set(gca,'FontSize',10);
xlim([100*min(FSDShapeGS_RelErr) 100*max(FSDShapeGS_RelErr)]);
ylim([100*min(meanFSDvar_GS) 100*max(meanFSDvar_GS)]);
title(sprintf('%s FSD Ground State Sampling (%u Samples) \nnormalization uncertainty %u %%',Isotopologue,nTrials, FSDNorm_RelErr*100));
grid on;
linkaxes([s1 s2 s3],'x');

save_dir = sprintf('./plots/FSD_diagnostic/');
save_name = sprintf('%s_Variance_Overview_GroundStates_%.0fSamples',Isotopologue,nTrials);
export_fig([save_dir,'png/',save_name,'.png']);
export_fig([save_dir,'fig/',save_name,'.fig']);
publish_figurePDF(fig68,[save_dir,'pdf/',save_name,'.pdf']);
%close;


