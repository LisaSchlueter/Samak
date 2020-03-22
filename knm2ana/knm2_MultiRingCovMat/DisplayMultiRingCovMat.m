% Script to Develop, Test and Display 
% Covariance matrices for multi-ring fit
% taking ring to ring correlation into account
% Lisa, Oct 2019

%% Set up multi ring model
RunList = 'KNM2_Prompt';
chi2 = 'chi2CMShape';
RingMerge = 'Full'; % 'Full == 4 rings
range =40;
E0 = knm2FS_GetE0Twins('SanityPlot','OFF');
% read data and set up model
CommonArg = {'RunList',RunList,...
    'chi2','chi2Stat',...
    'DataType','Twin',...
    'fixPar','mNu E0 Bkg Norm',...
    'RadiativeFlag','ON',...
    'minuitOpt','min ; migrad',...
    'SysBudget',33,...
    'AnaFlag','Ring',...
    'RingMerge',RingMerge,...
    'chi2',chi2,...
    'TwinBias_Q',E0,...
    'NonPoissonScaleFactor',1.112};

R = MultiRunAnalysis(CommonArg{:});
R.exclDataStart = R.GetexclDataStart(range);
%% Compute multi ring covariance matrix
R.ComputeCM;%
%('SysEffects',struct('TASR','ON','TCoff_OTHER','OFF','FSD','ON','RF_EL','ON','RF_BF','ON','RF_RX','ON','Stack','ON'),...
 %   'BkgCM','ON','nTrials',1000,'DataDriven','ON','RecomputeFlag','OFF')

%%
nrings = R.nRings;%numel(obj.SO.MACE_Ba_T); % always gives correct number of pseudo rings
nqUend = 25;
nqU_used = size(R.ModelObj.TBDIS(R.exclDataStart:nqUend,:),1);             % number of subruns, which are NOT excluded
exclIndex = sort(reshape(repmat(R.exclDataStart:nqUend,[nrings,1])+[0:nrings-1]'.*R.ModelObj.nqU,[nqU_used*nrings,1]));
close all
figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
[cp coeff]= corplot(R.FitCMFracShape(exclIndex,exclIndex));
pbaspect([1 1 1])
PrettyFigureFormat('FontSize',22);
cb = colorbar;
cb.Label.String = sprintf('Correlation coefficient \\rho');
cb.FontSize = 22;
set(gca,'XTick',[nqU_used/2+(0:nrings-1)*nqU_used]); set(gca,'YTick',[nqU_used/2+(0:nrings-1)*nqU_used]);
xlabel('Ring'); ylabel('Ring');
export_fig(gcf,'./plots/MultiRingCorrMat.pdf','-painters');

%%
close all
figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
imagesc(R.FitCMFrac(exclIndex,exclIndex));
pbaspect([1 1 1])
set(gca,'XTick',nqU_used/2+(0:nrings-1)*nqU_used); set(gca,'YTick',nqU_used/2+(0:nrings-1)*nqU_used);
set(gca,'XTickLabel',string(1:nrings))
set(gca,'YTickLabel',string(1:nrings))
PrettyFigureFormat('FontSize',24);
set(gca, 'XAxisLocation', 'top');
xlabel('Ring'); ylabel('Ring');
cb = colorbar;
cb.Label.String = sprintf('Fractional covariance');
cb.FontSize = 22;
export_fig(gcf,'./plots/MultiRingCovMatFrac.pdf','-painters');

