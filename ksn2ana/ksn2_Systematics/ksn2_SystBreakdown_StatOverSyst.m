% ksn2 extract systematic contribution
% raster scan
% compare variances
% systematics one-by-one
%% settings that might change
RasterScan = 'ON';
DataType = 'Twin';
nGridSteps = 30;
range = 40;
savedir = [getenv('SamakPath'),'ksn2ana/ksn2_Systematics/results/'];
savename = sprintf('%sksn2_SystBreakdown_StatOverSyst_%s_%.0feV_RasterScan%s.mat',savedir,DataType,range,RasterScan);
if exist(savename,'file')  
    load(savename);
     fprintf('load file %s \n',savename);
else
    
    %% configure RunAnalysis object
    RunAnaArg = {'RunList','KNM2_Prompt',...
        'DataType',DataType,...
        'fixPar','E0 Norm Bkg',...%free par
        'SysBudget',40,...
        'fitter','minuit',...
        'minuitOpt','min;migrad',...
        'RadiativeFlag','ON',...
        'FSDFlag','KNM2_0p5eV',...
        'ELossFlag','KatrinT2A20',...
        'AnaFlag','StackPixel',...
        'chi2','chi2Stat',...
        'NonPoissonScaleFactor',1,...
        'FSD_Sigma',sqrt(0.0124+0.0025),...
        'TwinBias_FSDSigma',sqrt(0.0124+0.0025),...
        'TwinBias_Q',18573.7,...
        'PullFlag',99,...;%99 = no pull
        'BKG_PtSlope',3*1e-06,...
        'TwinBias_BKG_PtSlope',3*1e-06,...
        'DopplerEffectFlag','FSD'};
    A = MultiRunAnalysis(RunAnaArg{:});
    %% configure Sterile analysis object
    SterileArg = {'RunAnaObj',A,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
        'nGridSteps',nGridSteps,...
        'SmartGrid','OFF',...
        'RecomputeFlag','OFF',...
        'SysEffect','all',...
        'RandMC','OFF',...
        'range',range,...
        'LoadGridArg',{'mNu4SqTestGrid',5}};
    
    %%
    S = SterileAnalysis(SterileArg{:});
    
    %%
    SysEffectsAll   = {'all','LongPlasma','BkgPT','NP','RF_BF','FSD', 'Bkg', 'RF_RX','FPDeff',...
        'TASR','RF_EL' ,'Stack','TCoff_OTHER'};
    SysEffectLabel    = {'Total'; 'Plasma';'Penning background';'Non-Poisson background';...
        'Magnetic fields';'Final-state distribution';'Background qU-slope';...
        'Source scattering';'Detector efficiency';
        'Tritium activity fluctuations';'Energy-loss function';...
        'HV fluctuations';'Theoretical corrections';};
    nSys = numel(SysEffectsAll);
    Ratio       = zeros(nSys,1e4);
    mNu4Sq      = zeros(nSys,1e4);
    sin2t4_Stat = zeros(nSys,1e4);
    sin2t4_Sys  = zeros(nSys,1e4);
    sin2t4_Tot  = zeros(nSys,1e4);
    
    for i=1:nSys
        S.SysEffect = SysEffectsAll{i};
        if strcmp(SysEffectsAll{i},'all') ||  strcmp(SysEffectsAll{i},'NP')
            S.RunAnaObj.NonPoissonScaleFactor = 1.112;
        else
            S.RunAnaObj.NonPoissonScaleFactor = 1;
        end
        
        [Ratio(i,:),~,mNu4Sq(i,:),sin2t4_Stat(i,:),sin2t4_Tot(i,:),sin2t4_Sys(i,:)] = S.StatOverSysKsn2('RasterScan',RasterScan);
        
    end
    % save
    MakeDir(savedir);
    save(savename,'Ratio','mNu4Sq','sin2t4_Stat','sin2t4_Tot',...
        'sin2t4_Sys','SysEffectsAll','SysEffectLabel','nSys',...
        'RunAnaArg','SterileArg'); 
    fprintf('save file to %s \n',savename);
end

%% plot 1. variance ratio
pHandle = cell(numel(SysEffectsAll),1);
LineStyle = {'-','-.',':','--','-','-.',':','--','-','-.',':','--','-.','-','--','-.',':','--','-.',':','--'};
Colors = colormap('jet');

fRatio = figure('Units','normalized','Position',[0.1,0.1,0.7,0.5]);
hold on;
for i=1:nSys
    if strcmp(SysEffectsAll{i},'all')
        PltColor = rgb('Black');
        LineWidth = 3;
    else
        PltColor =Colors(floor((i)*256/(nSys)),:);
        LineWidth = 2;
    end
    pHandle{i} = plot(sqrt(mNu4Sq(i,:)),Ratio(i,:),'LineWidth',LineWidth,'Color',PltColor,'LineStyle',LineStyle{i});
end
ylabel(sprintf('\\sigma_{syst.}^2 / \\sigma_{total}^2'))
xlabel(sprintf('{\\itm}_4^2 (eV^2)'));
set(gca,'XScale','log');
PrettyFigureFormat('FontSize',22);
xlim([1 40^2]);
 ylim([0 0.5])        
leg = legend([pHandle{:}],SysEffectLabel);
PrettyLegendFormat(leg);
leg.Title.String = 'Stystematic effect';
leg.Title.FontWeight ='normal';
leg.Location = 'eastoutside';


plotdir = strrep(savedir,'results','plots');
plotname = sprintf('%sksn2_SysBreakdown_StatOverSyst_Ratio_RasterScan%s.png',plotdir,RasterScan);
print(fRatio,plotname,'-dpng','-r400');
fprintf('save plot to %s \n',plotname);
export_fig(fRatio,strrep(plotname,'png','pdf'));
%% plot 2. only systematics contribution in contour plot style
pHandle = cell(numel(SysEffectsAll),1);
LineStyle = {'-','-.',':','--','-','-.',':','--','-','-.',':','--','-.','-','--','-.',':','--','-.',':','--'};
Colors = colormap('jet');

fSyst = figure('Units','normalized','Position',[0.1,0.1,0.7,0.5]);
hold on;
for i=1:nSys
    if strcmp(SysEffectsAll{i},'all')
        PltColor = rgb('Black');
        LineWidth = 3;
    else
        PltColor =Colors(floor((i)*256/(nSys)),:);
        LineWidth = 2;
    end
    pHandle{i} = plot(sin2t4_Sys(i,:),mNu4Sq(i,:),'LineWidth',LineWidth,'Color',PltColor,'LineStyle',LineStyle{i});
    
end
pStat = plot(sin2t4_Stat(1,:),mNu4Sq(1,:),'LineWidth',3,'Color',rgb('Silver'));
xlabel(sprintf('\\sigma_{syst.}'))
ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
set(gca,'XScale','log');
set(gca,'YScale','log');
PrettyFigureFormat('FontSize',22);
xlim([2e-04 0.5]);
         
leg = legend([pHandle{:},pStat],{SysEffectLabel{:},'Stat. only'});
PrettyLegendFormat(leg);
leg.Title.String = 'Stystematic effect';
leg.Title.FontWeight ='normal';
leg.Location = 'eastoutside';

plotdir = strrep(savedir,'results','plots');
plotname = sprintf('%sksn2_SysBreakdown_SystOnly_RasterScan%s.png',plotdir,RasterScan);
print(fSyst,plotname,'-dpng','-r400');
fprintf('save plot to %s \n',plotname);
export_fig(fSyst,strrep(plotname,'png','pdf'));