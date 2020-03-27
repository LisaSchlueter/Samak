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
    'minuitOpt','min ; minos',...
    'SysBudget',33,...
    'AnaFlag','Ring',...
    'RingMerge',RingMerge,...
    'chi2',chi2,...
    'TwinBias_Q',E0,...
    'NonPoissonScaleFactor',1.112};

R = MultiRunAnalysis(CommonArg{:});
R.exclDataStart = R.GetexclDataStart(range);
%% Compute multi ring covariance matrix
SysEffect = 'Bkg';
if strcmp(SysEffect,'all')
    R.ComputeCM;
elseif strcmp(SysEffect,'Bkg')
    R.ComputeCM('SysEffect',struct('FSD','OFF'))
else
    R.ComputeCM('SysEffect',struct(SysEffect,'ON'),'BkgCM','OFF')
end
%('SysEffects',struct('TASR','ON','TCoff_OTHER','OFF','FSD','ON','RF_EL','ON','RF_BF','ON','RF_RX','ON','Stack','ON'),...
%   'BkgCM','ON','nTrials',1000,'DataDriven','ON','RecomputeFlag','OFF')
%%
PlotStat = 'OFF'; %with or without stat. uncertainties
if strcmp(PlotStat,'OFF')
     CM          = R.FitCM_Obj.CovMat;
     CMFrac      = R.FitCM_Obj.CovMatFrac;
     CMFracShape = R.FitCM_Obj.CovMatFracShape;
else
     CM = R.FitCM;
     CMFrac = R.FitCMFrac;
     CMFracShape = R.FitCMFracShape;
end
%% correlation plot
savedir = [getenv('SamakPath'),'knm2ana/knm2_systematics/plots/'];
nrings = R.nRings;
nqUend = 32;
nqU_used = size(R.ModelObj.TBDIS(R.exclDataStart:nqUend,:),1);             % number of subruns, which are NOT excluded
exclIndex = sort(reshape(repmat(R.exclDataStart:nqUend,[nrings,1])+[0:nrings-1]'.*R.ModelObj.nqU,[nqU_used*nrings,1]));
close all
figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
[cp coeff]= corplot(CMFracShape(exclIndex,exclIndex));
pbaspect([1 1 1])
PrettyFigureFormat('FontSize',22);
cb = colorbar;
cb.Label.String = sprintf('Correlation coefficient \\rho');
cb.FontSize = 22;
set(gca,'XTick',[nqU_used/2+(0:nrings-1)*nqU_used]); set(gca,'YTick',[nqU_used/2+(0:nrings-1)*nqU_used]);
xlabel('Ring'); ylabel('Ring');
colormap(parula);
set(gca,'XMinorTick','off')
set(gca,'YMinorTick','off')
savename = sprintf('%sMultiRingCorrMatFracShape_%s_Stat%s.pdf',savedir,SysEffect,PlotStat);
export_fig(gcf,savename);
fprintf('save plot to %s \n',savename)
%%
close all
nqUend = 28;
nqU_used = size(R.ModelObj.TBDIS(R.exclDataStart:nqUend,:),1);             % number of subruns, which are NOT excluded
exclIndex = sort(reshape(repmat(R.exclDataStart:nqUend,[nrings,1])+[0:nrings-1]'.*R.ModelObj.nqU,[nqU_used*nrings,1]));

figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
imagesc(CMFracShape(exclIndex,exclIndex));
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
savename = sprintf('%sMultiRingCovMatFracShape_%s_Stat%s.pdf',savedir,SysEffect,PlotStat);
export_fig(gcf,savename);
fprintf('save plot to %s \n',savename)

%%  only 1 ring
plotdir = [getenv('SamakPath'),'knm2ana/knm2_systematics/plots/'];
R.FitCM_Obj.PlotCM('qUWindowIndexMax',200,'saveplot','ON','Convergence',...
    'OFF','Mode','Shape','PlotEffect','BkgSlopeRing1','savedir',plotdir,...
    'CovMatInput',CMFracShape(1:38,1:38));