

savedir = [getenv('SamakPath'),'knm1ana/knm1_systematics/knm1_covmats/'];%inputs/CovMat/Background/CM/'];
filename = sprintf('%sBackgroundCM_KNM1_293mcps_Constrain1e+09CpsPerEv_-5eVBkgRange_10000Trials.mat',savedir);
d = importdata(filename);
%%
qUStartIndex = 13;
MaxSlopeCpsPereV_1 = 15*1e-06;  %cut-off for linear fit (cov-mat)
InclLog1 = abs(d.Slopes)<=MaxSlopeCpsPereV_1;

MaxSlopeCpsPereV_2 = 5*1e-06;  %cut-off for linear fit (cov-mat)
InclLog2 = abs(d.Slopes)<=MaxSlopeCpsPereV_2; 
qU = d.obj.StudyObject.qU(qUStartIndex:end,1);
%% plot 1: randomized background spectra and fit slopes
close all
f1 = figure('Units','normalized','Position',[  0.0181    0.4044    0.6083    0.4700]);

Q_ref = 18574;

[pfit,afit] = boundedline(d.obj.StudyObject.qU(qUStartIndex:end,1)-Q_ref,...
    1e3.*nanmean(d.Bkg_Fit(qUStartIndex:end,:),2),...
    1e3.*std(d.Bkg_Fit(qUStartIndex:end,:),0,2,'omitnan'));
afit.FaceColor = rgb('Silver');
pfit.Color = rgb('DarkSlateGray');

hold on;

[pfit1,afit1] = boundedline(d.obj.StudyObject.qU(qUStartIndex:end,1)-Q_ref,...
    1e3.*nanmean(d.Bkg_Fit(qUStartIndex:end,InclLog1),2),...
    1e3.*std(d.Bkg_Fit(qUStartIndex:end,InclLog1),0,2,'omitnan'));
afit1.FaceColor = rgb('DodgerBlue');
pfit1.LineStyle = 'none';

[pfit2,afit2] = boundedline(d.obj.StudyObject.qU(qUStartIndex:end,1)-Q_ref,...
    1e3.*nanmean(d.Bkg_Fit(qUStartIndex:end,InclLog2),2),...
    1e3.*std(d.Bkg_Fit(qUStartIndex:end,InclLog2),0,2,'omitnan'));
afit2.FaceColor = rgb('SkyBlue');
pfit2.LineStyle = 'none';

for i=1:size(d.Data,1)
    sg = dscatter(squeeze(d.Data(i,1,:))-Q_ref+1e-09.*randn(1e4,1),1e3.*squeeze(d.Data(i,2,:)),'Msize',20); % density scatter plot
    hold on;
end

colormap(autumn);

set(gca,'FontSize',20);
xlabel(sprintf('Retarding energy qU - %.0f (eV)',Q_ref));
ylabel('Background (mcps)');
hold off;

leg = legend([afit,afit1,afit2],...% 'Randomized background spectra',...
    sprintf('All background slopes'),...
    sprintf('{\\its} (qU) \\leq %.0f mcps/keV',1e6*MaxSlopeCpsPereV_1),...
    sprintf('{\\its} (qU) \\leq   %.0f mcps/keV',1e6*MaxSlopeCpsPereV_2),...
  'Location','southwest');
leg.Title.String = sprintf('Linear fits 1\\sigma band');
leg.Title.FontWeight = 'normal';
PrettyLegendFormat(leg);
leg.FontSize = get(gca,'FontSize');

xlim([-41,d.obj.StudyObject.qU(end,1)-Q_ref+2]);
ylim([min(min(d.Data(:,2,:)))*1e3-0.005,max(max(d.Data(:,2,:)))*1e3+0.1])

PrettyFigureFormat('FontSize',22);

pltdir = [getenv('SamakPath'),'knm1ana/knm1_systematics/plots/'];
pltfile  =  [pltdir,sprintf('knm1_Bkg_syst_qU.pdf')];

export_fig(gcf,pltfile);
fprintf('Save plot to %s \n',pltfile);

%% plot 2: histogram fig2        = figure('Units','normalized','pos',[0.1 0.1 0.4 0.4]);
f1 = figure('Units','normalized','Position',[  0.0181    0.4044    0.4    0.4700]);

hall = histogram(d.Slopes.*1e6,'FaceColor',rgb('Silver'));
hold on;
SlopeOK = d.Slopes(logical(InclLog1));
h1 = histogram(d.Slopes(logical(InclLog1)).*1e06,'FaceColor',rgb('DodgerBlue'),...
    'FaceAlpha',1,'BinWidth',hall.BinWidth);
h2 = histogram(d.Slopes(logical(InclLog2)).*1e06,'FaceColor',rgb('SkyBlue'),...
    'FaceAlpha',1,'BinWidth',hall.BinWidth);
xlabel('Background slope (mcps / keV)');
ylabel('Occurrence');
PrettyFigureFormat('FontSize',22)

leg = legend([hall,h1,h2],...% 'Randomized background spectra',...
    sprintf('All background slopes'),...
    sprintf('{\\its} (qU) \\leq %.0f mcps/keV',1e6*MaxSlopeCpsPereV_1),...
    sprintf('{\\its} (qU) \\leq   %.0f mcps/keV',1e6*MaxSlopeCpsPereV_2),...
  'Location','northwest');
PrettyLegendFormat(leg);
leg.FontSize = get(gca,'FontSize');
ylim([0 750])
pltfile2  =  [pltdir,sprintf('knm1_Bkg_syst_qU_Hist.pdf')];

export_fig(gcf,pltfile2);
fprintf('Save plot to %s \n',pltfile2);

%% plot 3: covariance and correlation matrix
TBDIS = d.TBDIS_V;
TBDIS1 = d.TBDIS_V(:,logical(InclLog1));
TBDIS2 = d.TBDIS_V(:,logical(InclLog2));

%Compute Covariance Matrix
CovMat    = cov(TBDIS(qUStartIndex:end,:)');
CovMat1         = cov(TBDIS1(qUStartIndex:end,:)');
CovMat2         = cov(TBDIS2(qUStartIndex:end,:)');   

CovMatFrac = Cov2Frac(CovMat,d);
CovMatFrac1 = Cov2Frac(CovMat1,d);
CovMatFrac2 = Cov2Frac(CovMat2,d);
%%
close all
f1 = figure('Units','normalized','Position',[0.12 0.12,0.8,0.5]);
subplot(1,2,1);
imagesc(CovMatFrac1);
PrettyCovMat(qU);
 c = colorbar;
colormap(parula)
c.Label.String = sprintf('Fractional covariance');
c.Label.FontSize = get(gca,'FontSize')+4;
c.LineWidth = 1.5;
colormap(parula)
ax1 = gca;

subplot(1,2,2);
imagesc(corr(TBDIS1(qUStartIndex:end,:)'));
PrettyCovMat(qU);
c = colorbar;
c.Label.String = sprintf('Correlation coefficient');
c.Label.FontSize = get(gca,'FontSize')+4;
c.LineWidth = 1.5;
ax2 = gca;

ax1.Position(1) = 0.08;

pltfile3  =  [pltdir,sprintf('knm1_Bkg_syst_covmat.pdf')];

export_fig(gcf,pltfile3);
fprintf('Save plot to %s \n',pltfile3);


function cmfrac = Cov2Frac(cm,d)
  % Compute Fractional Covariance Matrix
BKGIS          = d.BKG_i.*d.obj.StudyObject.TimeSec.*d.obj.StudyObject.qUfrac(13:end);
cmfrac = cm./BKGIS./BKGIS';
cmfrac(isinf(cmfrac)) = 0; %for background region
cmfrac(isnan(cmfrac)) = 0; %for background region

end

function PrettyCovMat(q)
pbaspect([1 1 1])
ax = gca;
ax.XTick = 1:1:size(q,1);
ax.YTick = 1:1:size(q,1);


% xy ticks
TickLabels = strings(23,1);
for i=1:numel(q)
    if ismember(i,[1:7:20,21,27])
        TickLabels{i} = sprintf('%.0f eV',q(i)-18574);
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