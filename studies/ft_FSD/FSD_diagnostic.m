% Study to test FSD systematics treatment after Hamish Robertson's suggestion
% Introduce bin-to-bin uncorrelated fluctuation such that the variance of
% the FSD will have a variation (std) of 2-3%
Isotopologue = 'TT';
PlotFontSize = 18;
% Configuration
FSDNorm_RelErr = 0.01;%0.01;
FSDShapeGS_RelErr = [0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04 0.04];
FSDShapeES_RelErr = [0.04 0.05 0.07 0.1 0.15  0.2 0.25 0.3];%[0.01, 0.05, 0.1,  0.15, 0.2, 0.25];
nTrials = 10000;
state = 'excited state';

save_name = sprintf('./results/%s_FSDdiagnostic_%.0fSamples.mat',Isotopologue,nTrials);
if exist(save_name,'file')
    load(save_name);
else
% Set up Model 
A = MultiRunAnalysis('RunList','StackCD100all');

% Init Variables
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

FSD_P_norm = zeros(nProb,nTrials,numel(FSDShapeES_RelErr)); % Final State Probabilities after norm (happends inside ComputeTBDDS)
FSD_P      = zeros(nProb,nTrials,numel(FSDShapeES_RelErr)); % Final State Probabilities before norm
NormGS = zeros(nTrials,numel(FSDShapeES_RelErr));                            % Normalization Factor Ground State
NormES = zeros(nTrials,numel(FSDShapeES_RelErr));                            % Normalization Factor Excited State
NormBias = zeros(nTrials,numel(FSDShapeES_RelErr));

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


for ii=1:numel(FSDShapeES_RelErr)
% Do  Variation: Normalization
NormBias(:,ii) = randn(nTrials,1).*FSD_Norm_i*FSDNorm_RelErr;
% Do Variation: Bin-to-Bin uncorrelated
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

%variance distribution of samples widht different sampling uncertainties

FSDmean = zeros(nTrials,numel(FSDShapeES_RelErr));  % weighted mean energy   (eV)
FSDvar  = zeros(nTrials,numel(FSDShapeES_RelErr));  % weighted mean variance of FSD (eV^2)
FSDvar_GS  = zeros(nTrials,numel(FSDShapeES_RelErr));  % weighted mean variance of FSD (eV^2) for ground state
meanFSDvar = zeros(numel(FSDShapeES_RelErr),1);     % mean variance of FSD variance(=width)
meanFSDvar_GS = zeros(numel(FSDShapeES_RelErr),1);  
varFSDvar  = zeros(numel(FSDShapeES_RelErr),1);     % variance of FSD variance(=width)
varFSDvar_GS = zeros(numel(FSDShapeES_RelErr),1);

for ii=1:numel(FSDShapeES_RelErr)
for i=1:nTrials
    FSDmean(i,ii) = wmean(Energy,FSD_P_norm(:,i,ii)');
    FSDvar(i,ii) = var(Energy,FSD_P_norm(:,i,ii)');
    FSDvar_GS(i,ii) = var(Energy(1:GSESLimit), FSD_P_norm(1:GSESLimit,i,ii)');
end
meanFSDvar(ii)   = mean(FSDvar(:,ii));
meanFSDvar_GS(ii) = mean(FSDvar_GS(:,ii));
varFSDvar(ii)    = var(FSDvar(:,ii));
varFSDvar_GS(ii) = var(FSDvar_GS(:,ii));
end

FSD_P_norm_std  = squeeze(std(permute(FSD_P_norm,[2,1,3]))); % standard deviation of samles for plot
FSD_P_norm_mean = squeeze(mean(permute(FSD_P_norm,[2,1,3])));% mean of samles for plot
save(save_name)
end
%% FSD distribution for 1 sampling uncertainties
%for ii=1:numel(FSDShapeES_RelErr)
    IndexFSDRelErr=ii;
    IndexFSDRelErr=8;
mytitle = sprintf('%u DT Sample Distributions with \nnormalization uncertainty %u %% \n bin-to-bin uncorrelated uncertainties (ground state) %.1f \n different bin-to-bin uncorrelated uncertainties excited state',nTrials, FSDNorm_RelErr*100,FSDShapeGS_RelErr);
fig55 = figure('NumberTitle', 'off','Name',mytitle,'Renderer','opengl');
set(fig55, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.7]);
[l, p] = boundedline(Energy,FSD_P_theory*100,FSD_P_norm_std(:,IndexFSDRelErr)*100);%plot(Energy, FSD_P_norm(:,:,ii)*100,'--','Linewidth',1.5);
%l = plot(Energy,FSD_P_theory*100);
l.LineWidth = 2; l.Color = rgb('CadetBlue');
p.FaceColor = rgb('DarkOrange'); p.FaceAlpha = 0.5;
set(gca,'XScale','log')
xlabel('excitation energy (eV)');
ylabel('probability (%)');
title(sprintf('bin-to-bin uncorrelated %.1f%% (GS), %.1f%% (ES) \nnormalization %.1f%%', FSDShapeGS_RelErr(1)*100,FSDShapeES_RelErr(IndexFSDRelErr)*100,FSDNorm_RelErr*100'));
sample_leg = sprintf('1\\sigma error band');% );
theo_leg = sprintf('%s FSD Theory (%s)',Isotopologue, FSDCalc);
leg = legend([l p],theo_leg, sample_leg,'Location','northeast');
%leg = legend(l,theo_leg,'Location','northeast');
legend boxoff
PrettyFigureFormat;
xlim([min(Energy) max(Energy)]);
ylim([min(FSD_P_theory*100) 7.2]);%1.1*max(FSD_P_theory*100)]);
set(gca,'FontSize',PlotFontSize)
xlim([0.4 150]);

save_dir = sprintf('./plots/FSD_diagnostic/');
save_name = sprintf('%s_SampleFSD_AllStates_%.2gESerr_%.0fSamples',Isotopologue,FSDShapeES_RelErr(IndexFSDRelErr),nTrials);
print(gcf,[save_dir,'png/',save_name,'.png'],'-dpng','-r350');
export_fig([save_dir,'fig/',save_name,'.fig']);
publish_figurePDF(fig55,[save_dir,'pdf/',save_name,'.pdf'])
%close
%end

%% Delta FSD dstribution
for ii=1:numel(FSDShapeES_RelErr)
    IndexFSDRelErr=ii;
fig552 = figure('NumberTitle', 'off','Name',mytitle,'Renderer','opengl');
set(fig552, 'Units', 'normalized', 'Position', [0.9, 0.9, 0.8, 1.2]);
[l,p]=  boundedline(Energy,(FSD_P_theory-FSD_P_norm_mean(:,IndexFSDRelErr))*100,FSD_P_norm_std(:,IndexFSDRelErr)*100);
l.LineWidth = 1.5; l.Color = rgb('CadetBlue');
p.FaceColor = rgb('CadetBlue'); p.FaceAlpha = 0.5;
set(gca,'XScale','log')
xlabel('excitation energy (eV)');
ylabel('\Delta probability(%)');
sample_title = sprintf('bin-to-bin uncorrelated uncertainty %.1f%% (GS), %.1f%% (ES)  \nnormalization uncertainty %.1f%%', FSDShapeGS_RelErr(1)*100,FSDShapeES_RelErr(IndexFSDRelErr)*100,FSDNorm_RelErr*100);
title(sample_title);
sample_leg = sprintf('1\\sigma error band');% );
theo_leg = sprintf('%s Theory FSD (%s)',Isotopologue, FSDCalc);
leg = legend([l p],theo_leg, sample_leg,'Location','northwest');
legend boxoff
PrettyFigureFormat;
xlim([min(Energy) max(Energy)]);
ylim([-2, 1.1])
%ylim([min((FSD_P_theory-FSD_P_norm_mean(:,IndexFSDRelErr))*100)*1.2-max(FSD_P_norm_std(:,IndexFSDRelErr))*100, max(FSD_P_norm_std(:,IndexFSDRelErr))*100+max((FSD_P_theory-FSD_P_norm_mean(:,IndexFSDRelErr))*100)])
set(gca,'FontSize',PlotFontSize)


save_dir = sprintf('./plots/FSD_diagnostic/');
save_name = sprintf('%s_DeltaSampleFSD_AllStates_%.2gESerr_%.0fSamples',Isotopologue,FSDShapeES_RelErr(IndexFSDRelErr),nTrials);
export_fig([save_dir,'png/',save_name,'.png']);
export_fig([save_dir,'fig/',save_name,'.fig']);
publish_figurePDF(fig552,[save_dir,'pdf/',save_name,'.pdf']);
 close;
end
%% Variance Distribution
fig67 = figure('NumberTitle', 'off','Name','variance','Renderer','opengl');
set(fig67, 'Units', 'normalized', 'Position', [0.9, 0.9, 0.7, 0.7]);
IndexFSDRelErr= 6;
h1 = nhist(FSDvar(:,IndexFSDRelErr),'color',rgb('CadetBlue'));
sample_leg = sprintf('\nstd/mean=%.1f %%',sqrt(varFSDvar(IndexFSDRelErr))./meanFSDvar(IndexFSDRelErr)*100);
l1 = legend(sprintf('%s FSD',Isotopologue), [h1(10:31),sample_leg],'Location','northeast');%,['\sigma_{theo} = ',theo_leg],['<\sigma²>  = ',var_leg1,'var(\sigma²)= ',var_leg2]);
legend boxoff
xlabel('variance \sigma^2 (eV²)'); ylabel('');
set(gca,'FontSize',PlotFontSize);
PrettyFigureFormat;
l1.FontSize = PlotFontSize;
title(sprintf('bin-to-bin uncorrelated uncertainty GS %.1f%% - ES %.1f%% , normalization uncertainty %.1f%%',...
    FSDShapeGS_RelErr(1)*100,FSDShapeES_RelErr(IndexFSDRelErr)*100,FSDNorm_RelErr*100),'FontSize',PlotFontSize-2);
       
save_dir = sprintf('./plots/FSD_diagnostic/');
save_name = sprintf('%s_Variance_AllStates_Eserr%.2g_%.0fSamples',Isotopologue,FSDShapeES_RelErr(IndexFSDRelErr),nTrials);
export_fig([save_dir,'png/',save_name,'.png']);
export_fig([save_dir,'fig/',save_name,'.fig']);
publish_figurePDF(fig67,[save_dir,'pdf/',save_name,'.pdf']);
close;
%% std of width of FSD depending in sampling uncertainty
fig682 = figure('Name','FSDstd1','Renderer','opengl');
set(fig682, 'Units', 'normalized', 'Position', [0.9, 0.9, 0.7, 0.7]);

plot(FSDShapeES_RelErr*100,sqrt(varFSDvar)./meanFSDvar*100,'--o','MarkerSize',8,'LineWidth',3,'Color',rgb('FireBrick'),'MarkerFaceColor',rgb('FireBrick'));
hold on;
plot(FSDShapeES_RelErr*100,2*ones(numel(FSDShapeES_RelErr),1),'--k','LineWidth',2);
plot(FSDShapeES_RelErr*100,3*ones(numel(FSDShapeES_RelErr),1),'--k','LineWidth',2);
xlabel(sprintf('FSD bin-to-bin uncorrelated uncertainty (%s)',state));
ylabel(['\sigma ',sprintf('of FSD variance / mean variance (%%)')]);
PrettyFigureFormat;
set(gca,'FontSize',PlotFontSize);
grid on;
xlim([100*min(FSDShapeES_RelErr) 100*max(FSDShapeES_RelErr)]);
ylim([min(sqrt(varFSDvar)./meanFSDvar*100) max(sqrt(varFSDvar)./meanFSDvar*100)]);
title(sprintf('%s FSD Sampling (%u Samples) \nnormalization uncertainty %u %% + bin-to-bin uncorr. uncertainty %u %% (GS)',...
    Isotopologue, nTrials, FSDNorm_RelErr*100,FSDShapeGS_RelErr(1)*100),'FontSize',PlotFontSize-2);


save_dir = sprintf('./plots/FSD_diagnostic/');
save_name = sprintf('%s_Variance_ESErr_AllStates_%.0fSamples',Isotopologue,nTrials);
export_fig([save_dir,'png/',save_name,'.png']);
export_fig([save_dir,'fig/',save_name,'.fig']);
publish_figurePDF(fig682,[save_dir,'pdf/',save_name,'.pdf']);
%% overview std of width of FSD depending in sampling uncertainty
fig68 = figure('Name','FSDstd','Renderer','opengl');
set(fig68, 'Units', 'normalized', 'Position', [0.9, 0.9, 0.5, 1]);

s1 = subplot(3,1,1);
plot(FSDShapeES_RelErr*100,sqrt(varFSDvar),'--o','MarkerSize',8,'LineWidth',3,'Color',rgb('IndianRed'));
xlabel(sprintf('FSD bin-to-bin uncorrelated uncertainty (%s)',state));
ylabel('\sigma of FSD variance (eV^2)')
PrettyFigureFormat;
set(gca,'FontSize',12);
xlim([100*min(FSDShapeES_RelErr) 100*max(FSDShapeES_RelErr)]);
ylim([min(sqrt(varFSDvar)) max(sqrt(varFSDvar))]);
title(sprintf('%s FSD Sampling (%u Samples) \nnormalization uncertainty %u %% + bin-to-bin uncorr. uncertainty %u %% (GS)',Isotopologue, nTrials, FSDNorm_RelErr*100,FSDShapeGS_RelErr(1)*100));

s2 = subplot(3,1,2);
plot(FSDShapeES_RelErr*100,meanFSDvar,'--o','MarkerSize',8,'LineWidth',3,'Color',rgb('IndianRed'));
hold on;
xlabel(sprintf('FSD bin-to-bin uncorrelated uncertainty (%s)',state));
ylabel(sprintf('mean variance (eV^2)'));
PrettyFigureFormat;
set(gca,'FontSize',12);
xlim([100*min(FSDShapeES_RelErr) 100*max(FSDShapeES_RelErr)]);
ylim([min(meanFSDvar) max(meanFSDvar)]);
linkaxes([s1,s2],'x');

s2 = subplot(3,1,3);
plot(FSDShapeES_RelErr*100,sqrt(varFSDvar)./meanFSDvar*100,'--o','MarkerSize',8,'LineWidth',3,'Color',rgb('IndianRed'));
hold on;
plot(FSDShapeES_RelErr*100,2*ones(numel(FSDShapeES_RelErr),1),'--k');
plot(FSDShapeES_RelErr*100,3*ones(numel(FSDShapeES_RelErr),1),'--k');
xlabel(sprintf('FSD bin-to-bin uncorrelated uncertainty (%s)',state));
ylabel(['\sigma ',sprintf('/ mean variance (%%)')]);
PrettyFigureFormat;
set(gca,'FontSize',12);
xlim([100*min(FSDShapeES_RelErr) 100*max(FSDShapeES_RelErr)]);
ylim([min(sqrt(varFSDvar)./meanFSDvar*100) max(sqrt(varFSDvar)./meanFSDvar*100)]);
linkaxes([s1,s2],'x');

save_dir = sprintf('./plots/FSD_diagnostic/');
save_name = sprintf('%s_Variance_Overview_AllStates_%.0fSamples',Isotopologue,nTrials);
export_fig([save_dir,'png/',save_name,'.png']);
export_fig([save_dir,'fig/',save_name,'.fig']);
publish_figurePDF(fig68,[save_dir,'pdf/',save_name,'.pdf']);

%%
fig69 = figure('Name','FSDstd','Renderer','opengl');
set(fig69, 'Units', 'normalized', 'Position', [0.9, 0.9, 0.7, 0.7]);

d = importdata('../../inputs/CovMat/FSD/CM/FSD_DT-CovMat_1000Trials_StackCD100allex2_0.01NormErr_0.04GS_0.18ES_ShapeErr.mat');
plot(A.ModelObj.qU-18575,d.TBDIS_V,'LineWidth',3)
PrettyFigureFormat;
xlabel('retarding potential - E_0 (eV)');
ylabel('counts');
set(gca,'FontSize',24);
xlim([-200,20]);
ylim([0 5.4*1e5]);


%print(fig69,'./plots/FSD_CM_Spectrum','-dpng','-r350');
