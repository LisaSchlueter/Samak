% ksn2 systematics breakdown
% combi plot with syst. only and ratio
SavePlt = 'ON';
RasterScan = 'ON';
DataType = 'Twin';
nGridSteps = 30;
range = 40;
InterpMode = 'spline';
savedir = [getenv('SamakPath'),'ksn2ana/ksn2_Systematics/results/'];
savename = sprintf('%sksn2_SystBreakdown_StatOverSyst_%s_%.0feV_RasterScan%s_%sInterp.mat',...
    savedir,DataType,range,RasterScan,InterpMode);
if exist(savename,'file')
    d = importdata(savename);
    fprintf('load file %s \n',savename);
else
    fprintf('file not found %s \n Run ksn2_SystBreakdown_StatOverSyst.m',savename);
end

%%
FontSize = 24;
pHandle = cell(d.nSys,1);

f1 = figure('Units','normalized','Position',[-0.1,-0.1,0.7,0.6]);
s1 = subplot(1,4,[1:2]);
hold on;
for i=1:d.nSys
    if strcmp(d.SysEffectsAll{i},'all')
        PltColor = rgb('Black');
        LineWidth = 3;
    else
        PltColor =Colors(floor((i)*256/(d.nSys)),:);
        LineWidth = 2;
    end
    pHandle{i} = plot(d.sin2t4_Sys(i,:),d.mNu4Sq(i,:),'LineWidth',LineWidth,'Color',PltColor,'LineStyle',LineStyle{i});
    
end
pStat = plot(d.sin2t4_Stat(1,:),d.mNu4Sq(1,:),'LineWidth',3,'Color',rgb('Silver'));
xlabel(sprintf('\\sigma_{syst.}'))
ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
xlim([2e-04 0.5]);
ylim([1.3 38^2])
PrettyFigureFormat('FontSize',FontSize);
set(gca,'XScale','log');
set(gca,'YScale','log');
ax1 = gca;
xlim([3e-04 0.5])

% legend
SysEffectLabel = {'All systematics effects',...
               'Source-potential variations',...
           'Scan-step-duration-dependent background',...
           'Non-Poisson background',...
    'Magnetic fields',...
    'Molecular final-state distribution',...
    sprintf('{\\itqU}-dependent background'),...
    sprintf('Column density \\times inelastic scat. cross-section'),...
    'Detector efficiency',...
    'Activity fluctuations',...
    'Energy-loss function',...
    'High voltage stability and reproducibility',...
    'Theoretical corrections','Statistics only'};
leg = legend([pHandle{:},pStat],SysEffectLabel);
PrettyLegendFormat(leg);
% leg.Title.String = 'Stystematic effect';
% leg.Title.FontWeight ='normal';
leg.Location = 'northoutside';%'eastoutside';

%% ratio subplot 2
s2 = subplot(1,4,[3:4]);
plot(0.5.*ones(10,1),linspace(0.1,2e3,10),'k-','LineWidth',1.5)
t = text(0.45,40,'stat. = syst.','HorizontalAlignment','center','Rotation',90,'FontSize',18,'FontWeight','normal');
hold on
for i=1:nSys
    if strcmp(SysEffectsAll{i},'all')
        PltColor = rgb('Black');
        LineWidth = 3;
    else
        PltColor =Colors(floor((i)*256/(nSys)),:);
        LineWidth = 2;
    end
    plot(sin2t4_Sys(i,:).^2./sin2t4_Tot(1,:).^2,mNu4Sq(i,:),'LineWidth',LineWidth,'Color',PltColor,'LineStyle',LineStyle{i});
end
pStat = plot(sin2t4_Stat(1,:).^2./sin2t4_Tot(1,:).^2,mNu4Sq(1,:),...
    'LineWidth',3,'Color',rgb('Silver'),'LineStyle','-');
xlabel(sprintf('\\sigma^2 / \\sigma_{total}^2'))
PrettyFigureFormat('FontSize',FontSize);
set(gca,'XScale','lin');
set(gca,'YScale','log');
linkaxes([s1,s2],'y');
ax2 = gca;
ax2.YAxisLocation = 'right';
%yticklabels('');
xlim([0 1])
ylim([1 1600])
%% positions
ax2.Position(4) = 0.55;
ax1.Position(4) = 0.55;
ax1.Position(1) = 0.1;
ax1.Position(2) = 0.15;
ax2.Position(2) = 0.15;

ax1.Position(3) = 0.5;
ax2.Position(3) = 0.19;
ax2.Position(1) = 0.61;

leg.NumColumns = 2;
leg.Position(1) = 0.03;

plotdir = [];
