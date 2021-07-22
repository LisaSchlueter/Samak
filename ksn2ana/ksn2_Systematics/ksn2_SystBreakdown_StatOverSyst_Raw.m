% ksn2 extract systematic contribution
% raster scan
% compare variances
% systematics one-by-one
% without interpolation 
%% settings that might change
SavePlt = 'OFF';
RasterScan = 'ON';
DataType = 'Twin';
nGridSteps = 30;
range = 40;

savedir = [getenv('SamakPath'),'ksn2ana/ksn2_Systematics/results/'];
savename = sprintf('%sksn2_SystBreakdown_StatOverSyst_Raw_%s_%.0feV_RasterScan%s.mat',...
    savedir,DataType,range,RasterScan);
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
    if strcmp(DataType,'Real')
    LoadGridArg = {'mNu4SqTestGrid',5,'ExtmNu4Sq','ON'};
else
    LoadGridArg = {'mNu4SqTestGrid',5};
    end

    SterileArg = {'RunAnaObj',A,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
        'nGridSteps',nGridSteps,...
        'SmartGrid','OFF',...
        'RecomputeFlag','OFF',...
        'SysEffect','all',...
        'RandMC','OFF',...
        'range',range,...
        'LoadGridArg',LoadGridArg};
    
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

    % stat only
    S.RunAnaObj.chi2 = 'chi2Stat';
    S.RunAnaObj.NonPoissonScaleFactor = 1;
    S.LoadGridFile(S.LoadGridArg{:});
    S.chi2 = S.chi2';
    S.ContourPlot('RasterScan',RasterScan); close;
    mNu4Sq = S.mNu4Sq_contour;
    sin2t4_Stat = S.sin2T4_contour;
   
    sin2t4_Sys  = zeros(nSys,numel(sin2t4_Stat));
    sin2t4_Tot  = zeros(nSys,numel(sin2t4_Stat));
    
    for i=1:nSys
        S.SysEffect = SysEffectsAll{i};
        if strcmp(SysEffectsAll{i},'all') ||  strcmp(SysEffectsAll{i},'NP')
            S.RunAnaObj.NonPoissonScaleFactor = 1.112;
        else
            S.RunAnaObj.NonPoissonScaleFactor = 1;
        end
        
        if ~strcmp(S.SysEffect,'NP')
            S.RunAnaObj.chi2 = 'chi2CMShape';
        else
              S.RunAnaObj.chi2 = 'chi2Stat';
        end
        
        S.LoadGridFile(S.LoadGridArg{:});
        S.chi2 = S.chi2';
        S.ContourPlot('RasterScan',RasterScan); close;
        mNu4Sq_Tot = S.mNu4Sq_contour;
        sin2t4_Tot_tmp = S.sin2T4_contour;
        if i==nSys
            a =1
        end
        % match stat. only
        sin2t4_Tot(i,:) = interp1(mNu4Sq_Tot,sin2t4_Tot_tmp,mNu4Sq,'lin');
        sin2t4_Sys(i,:) = sqrt(sin2t4_Tot(i,:).^2-sin2t4_Stat.^2);
    end
    
   Ratio = sin2t4_Sys.^2./sin2t4_Tot.^2;
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
    pHandle{i} = plot(mNu4Sq,sin2t4_Sys(i,:).^2./sin2t4_Tot(1,:).^2,'LineWidth',LineWidth,'Color',PltColor,'LineStyle',LineStyle{i});
end
pStat = plot(mNu4Sq,sin2t4_Stat.^2./sin2t4_Tot(1,:).^2,'LineWidth',LineWidth,'Color',rgb('Silver'),'LineStyle','-');
ylabel(sprintf('\\sigma_{syst.}^2 / \\sigma_{total}^2'))
xlabel(sprintf('{\\itm}_4^2 (eV^2)'));
set(gca,'XScale','log');
PrettyFigureFormat('FontSize',22);
xlim([1.6 40^2]);
 ylim([0 0.5])        
leg = legend([pHandle{:}],SysEffectLabel);
PrettyLegendFormat(leg);
leg.Title.String = 'Systematic effect';
leg.Title.FontWeight ='normal';
leg.Location = 'eastoutside';

if strcmp(SavePlt,'ON')
    plotdir = strrep(savedir,'results','plots');
    plotname = sprintf('%sksn2_SysBreakdown_StatOverSyst_Ratio_RasterScan_%s_%s_%sInterp.png',...
        plotdir,DataType,RasterScan,InterpMode);
    print(fRatio,plotname,'-dpng','-r400');
    fprintf('save plot to %s \n',plotname);
    export_fig(fRatio,strrep(plotname,'png','pdf'));
end

%%
for i=1:nSys
fprintf('%.0f mean = %.3f \n',i,nanmean(Ratio(i,:)));
end