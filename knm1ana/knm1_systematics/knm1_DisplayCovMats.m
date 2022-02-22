close all
clear all
SysEffect = 'FPD';

savedir = [getenv('SamakPath'),'knm1ana/knm1_systematics/knm1_covmats/'];
switch SysEffect
    case 'NP'
        filename = sprintf('%sknm1_NonPoisCM.mat',savedir);
    case 'FSD'
        filename = sprintf('%sFSD_SibilleFull_DT-HT-TT-CovMat_5000Trials_KNM1_0.01NormErr_0.04GS_0.18ES_ShapeErr.mat',savedir);
    case 'TASR'
        filename = sprintf('%sTASR_CovMat_KNM1_7.281e-09RelErr.mat',savedir);
    case 'Bfield'
        filename = sprintf('%sWGTSMACE_CovMat_KNM1_CD1.10992e+17molPercm2_IsXEdep_Bs2.52T-0.025err_Bmax4.23T-0.002err_Ba6.00T-0.01err_KatrinT2-ELossOFF_elossBinning-9e+03eV-0.2eVStep_RFbinStep-0.04eV_1000Trials.mat',savedir);
    case 'RhoD'
        filename = sprintf('%sWGTSMACE_CovMat_KNM1_CD1.10992e+17molPercm2-0.0085err_IsXEdep-0err_Bs2.52T_Bmax4.23T_Ba6.00T_KatrinT2-ELossOFF_elossBinning-9e+03eV-0.2eVStep_RFbinStep-0.04eV_1000Trials.mat',savedir);
    case 'ELoss'
        filename = sprintf('%sWGTSMACE_CovMat_KNM1_CD1.10992e+17molPercm2_IsXEdep_Bs2.52T_Bmax4.23T_Ba6.00T_KatrinT2_elossBinning-5e+02eV-0.2eVStep_RFbinStep-0.04eV_1000Trials.mat',savedir);
    case 'TC'
        filename = sprintf('%sTCoff_CovMat_KNM1.mat',savedir);
    case 'Tot'
        filename = sprintf('%sknm1_TotalFitCM.mat',savedir); 
    case 'FPD'
       filename = sprintf('%sFPDeff_CovMat_KNM1_0.0001-RelErr_5000-Trials.mat',savedir);     
end

d = importdata(filename);

qUMax = 39-8;

if ~strcmp(SysEffect,'NP') && ~strcmp(SysEffect,'Tot')
    qU = d.obj.StudyObject.qU;%(13:qUMax);
    IdxBkg = find(qU-18574>0);
    PltCM   = d.obj.CovMatFracShape;
    CorrMat = corrcov(d.obj.CovMatFrac);
    PltCM(IdxBkg,IdxBkg) = 0;
    CorrMat(isnan(CorrMat)) = 0; %for background region
elseif strcmp(SysEffect,'NP')
    qU = d.qU;
    IdxBkg = find(qU-18574>0);
    PltCM   = d.NP_CMFrac;
    CorrMat = corrcov(d.NP_CMFrac);
    PltCM(IdxBkg,IdxBkg) = 0;
    CorrMat(isnan(CorrMat)) = 0; %for background region
elseif strcmp(SysEffect,'Tot')
    qU = d.qU;
    IdxBkg = find(qU-18574>0);
    PltCM   = d.FitCMFracShape;
    CorrMat = corrcov(d.FitCM);
    PltCM(IdxBkg,IdxBkg) = 0;
    CorrMat(isnan(CorrMat)) = 0; %for background region
end
%%
close all
f1 = figure('Units','normalized','Position',[0.12 0.12,0.8,0.5]);
subplot(1,2,1);
imagesc(PltCM(13:qUMax,13:qUMax));
PrettyCovMat(qU(13:qUMax));
c = colorbar;
colormap(parula)
c.Label.String = sprintf('Fractional covariance');
c.Label.FontSize = get(gca,'FontSize')+4;
c.LineWidth = 1.5;
colormap(parula)
ax1 = gca;

subplot(1,2,2);
imagesc(CorrMat(13:end,13:end));
PrettyCovMat(qU(13:end));
c = colorbar;
c.Label.String = sprintf('Correlation coefficient');
c.Label.FontSize = get(gca,'FontSize')+4;
c.LineWidth = 1.5;
ax2 = gca;
ax1.Position(1) = 0.08;

%%
pltdir   = [getenv('SamakPath'),'knm1ana/knm1_systematics/plots/'];
pltfile  =  sprintf('%sknm1_covmat_%s.pdf',pltdir,SysEffect);
export_fig(gcf,pltfile);
fprintf('Save plot to %s \n',pltfile);


function PrettyCovMat(q)
pbaspect([1 1 1])
ax = gca;
ax.XTick = 1:1:size(q,1);
ax.YTick = 1:1:size(q,1);


% xy ticks
TickLabels = strings(23,1);
for i=1:numel(q)
    if numel(q)<20
        if ismember(i,[1:6:20,numel(q)])
            TickLabels{i} = sprintf('%.0f eV',q(i)-18574);
        end
    else
        if ismember(i,[1:7:30,numel(q)])
            TickLabels{i} = sprintf('%.0f eV',q(i)-18574);
        end
    end
end

xlabel('Retarding energy - 18574 eV')
ylabel('Retarding energy - 18574 eV')
PrettyFigureFormat('FontSize',20);
set(gca,'XMinorTick','off');
set(gca,'YMinorTick','off');
xticklabels(TickLabels);
yticklabels(TickLabels);

end

