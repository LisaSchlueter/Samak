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
LocalFontSize = 16;
pHandle = cell(d.nSys,1);
LineStyle = {'-','-.',':','--','-','-.',':','--','-','-.',':','--','-.','-','--','-.',':','--','-.',':','--'};
Colors = colormap('jet');

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
xlabel(sprintf('|{\\itU}_{e4}|^2 sensitivity at 68.3%% C.L. '))%\\sigma_{syst.}
ylabel(sprintf('{\\itm}_4^2 (eV^{ 2})'));
xlim([2e-04 0.5]);
ylim([1.3 38^2])
PRLFormat;
set(gca,'FontSize',LocalFontSize);
set(gca,'XScale','log');
set(gca,'YScale','log');
ax1 = gca;
xlim([3e-04 0.5])

% legend
SysEffectLabel = {'Statistical uncertainty',...1
    'Combined systematic uncertainties',...2
    'Non-Poisson background',...3
    'Source-potential variations',...
    'Scan-step-duration-dependent background',...
    'Magnetic fields',...
    'Molecular final-state distribution',...
    sprintf('{\\itqU}-dependent background'),...
    sprintf('Column density \\times inelastic scat. cross-section'),...
    'High voltage stability and reproducibility',... % 12
    'Activity fluctuations',...
    'Energy-loss function',...
    'Detector efficiency',...% 9
    'Theoretical corrections'};
leg = legend([pStat,pHandle{1},pHandle{4},pHandle{[2,3,5:8]},pHandle{12},pHandle{[9,10,11,13]}],SysEffectLabel);
PrettyLegendFormat(leg);
% leg.Title.String = 'Stystematic effect';
% leg.Title.FontWeight ='normal';
leg.Location = 'northoutside';%'eastoutside';
leg.ItemTokenSize = [20 18];
%% ratio subplot 2
s2 = subplot(1,4,[3:4]);
plot(0.5.*ones(10,1),linspace(0.1,2e3,10),'k:','LineWidth',1.5)
t = text(0.45,40,'stat. = syst.','HorizontalAlignment','center','Rotation',90,'FontSize',LocalFontSize,'FontWeight','normal',...
    'FontName','Times New Roman');
hold on
for i=1:d.nSys
    if strcmp(d.SysEffectsAll{i},'all')
        PltColor = rgb('Black');
        LineWidth = 3;
    else
        PltColor =Colors(floor((i)*256/(d.nSys)),:);
        LineWidth = 2;
    end
    plot(d.sin2t4_Sys(i,:).^2./d.sin2t4_Tot(1,:).^2,d.mNu4Sq(i,:),'LineWidth',LineWidth,'Color',PltColor,'LineStyle',LineStyle{i});
end
pStat = plot(d.sin2t4_Stat(1,:).^2./d.sin2t4_Tot(1,:).^2,d.mNu4Sq(1,:),...
    'LineWidth',3,'Color',rgb('Silver'),'LineStyle','-');
xlabel(sprintf('\\sigma^2 / \\sigma_{total}^2'))
PRLFormat;
set(gca,'FontSize',LocalFontSize);
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
ax2.Position(3) = 0.22;%19;
ax2.Position(1) = 0.61;

leg.NumColumns = 3;
leg.Position(1) = 0.058;%3;
leg.FontSize = LocalFontSize; %16.5

ax1.XLabel.FontSize = LocalFontSize;
ax1.YLabel.FontSize = LocalFontSize;
ax2.XLabel.FontSize = LocalFontSize;
%%
plotdir = [getenv('SamakPath'),'ksn2ana/ksn2_Systematics/plots/'];
pltname = sprintf('%sksn2_SystBreakdown_%s_CombiPlt.pdf',plotdir,DataType);
export_fig(pltname);
%pltname = sprintf('%sksn2_SystBreakdown_%s_CombiPlt.png',plotdir,DataType);
%print(pltname,'-dpng','-r350');
