% ksn2 calculate chi2 grid search
% compare m2 free, m2 nuisance parameter
%% settings that might change
chi2 = 'chi2Stat';
DataType = 'Twin';
nGridSteps = 30;
range = 40;
%% configure RunAnalysis object
if strcmp(chi2,'chi2Stat')
    NonPoissonScaleFactor = 1;
elseif  strcmp(chi2,'chi2CMShape')
    NonPoissonScaleFactor = 1.112;
end

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
    'chi2',chi2,...
    'NonPoissonScaleFactor',NonPoissonScaleFactor,...
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
S.InterpMode = 'spline';
S.RunAnaObj.chi2 = 'chi2CMShape';
S.RunAnaObj.NonPoissonScaleFactor = 1;
SysEffectsAll   = {'Stat','LongPlasma','BkgPT','NP','FSD', 'RF_BF','Bkg','TASR','RF_EL',...
                     'RF_RX','FPDeff','Stack','TCoff_OTHER','all'}; 

HoldOn = 'OFF';
pHandle = cell(numel(SysEffectsAll),1);
LineStyle = {'-','-.',':','--','-','-.',':','--','-','-.',':','--','-.','-','--','-.',':','--','-.',':','--'};
Colors = colormap('jet');

Contour_x = cell(numel(SysEffectsAll),1);
Contour_y = cell(numel(SysEffectsAll),1);

m4Contour = logspace(log10(1.5),log10(34^2),1e3);
s4Contour = zeros(numel(SysEffectsAll),1e3);
for i=1:numel(SysEffectsAll)
    S.SysEffect = SysEffectsAll{i};
    LineWidth = 2;
    myColor = Colors(i*floor(256/numel(SysEffectsAll)),:);
    
    if strcmp(SysEffectsAll{i},'NP') || strcmp(SysEffectsAll{i},'Stat')
        S.RunAnaObj.chi2 = 'chi2Stat';
        if strcmp(SysEffectsAll{i},'NP')
            S.RunAnaObj.NonPoissonScaleFactor = 1.112;
        else
            LineWidth = 1.5;
            myColor = rgb('Silver');
        end
        
        S.LoadGridFile(S.LoadGridArg{:});
        S.Interp1Grid('Maxm4Sq',38^2);
        
        S.RunAnaObj.chi2 = 'chi2CMShape';
        S.RunAnaObj.NonPoissonScaleFactor = 1;
    elseif strcmp(SysEffectsAll{i},'all')
        S.RunAnaObj.NonPoissonScaleFactor = 1.112;
         S.LoadGridFile(S.LoadGridArg{:});
        S.Interp1Grid('Maxm4Sq',38^2);
        LineWidth = 2;
         myColor = rgb('Black');
    else
        
        S.LoadGridFile(S.LoadGridArg{:});
        S.Interp1Grid('Maxm4Sq',38^2);
    end

    pHandle{i} =  S.ContourPlot('HoldOn',HoldOn,'SavePlot','OFF','Color',myColor,...
        'LineStyle',LineStyle{i});
    pHandle{i}.LineWidth = LineWidth;
    HoldOn = 'ON';
    Contour_x{i} = S.sin2T4_contour;
    Contour_y{i} = S.mNu4Sq_contour;
    try
        s4Contour(i,:) = interp1(S.mNu4Sq_contour,S.sin2T4_contour,m4Contour,'spline');
    catch
        fprintf(2,'%s interp fail - try with reduced m4^2 \n',SysEffectsAll{i})
        S.sin2T4_contour = S.sin2T4_contour(6:end);
        S.mNu4Sq_contour = S.mNu4Sq_contour(6:end);
        s4Contour(i,:) = interp1(S.mNu4Sq_contour,S.sin2T4_contour,m4Contour,'spline');
    end
end

xlim([7e-03,0.5])
SysEffectLeg    = {'Stat. only'; 'Plasma';'Penning background';'Non-Poisson background';...
                     'Final-state distribution'; 'Magnetic fields';'Background qU-slope';...
                     'Tritium activity fluctuations';'Energy-loss function';...
                    'Source scattering';'Detector efficiency';...
                     'HV fluctuations';'Theoretical corrections';'Stat. and total syst.'};
               
leg = legend([pHandle{:}],SysEffectLeg);
PrettyLegendFormat(leg);
leg.Title.String = 'Stat. + 1 syst.';
leg.Title.FontWeight ='normal';
leg.Location = 'east';
% save 
plotname123 = sprintf('%sBudget%.0f_SystBreakdown.png',extractBefore(S.DefPlotName,'Budget'),A.SysBudget);
print(gcf,plotname123,'-dpng','-r300');
fprintf('save plot to %s \n',plotname123)
%% save zoom
ylim([100 800])
xlim([8.5e-03 0.017])
set(gca,'XScale','lin');
plotname123 = sprintf('%sBudget%.0f_SystBreakdownZoom.png',extractBefore(S.DefPlotName,'Budget'),A.SysBudget);
print(gcf,plotname123,'-dpng','-r300');
fprintf('save plot to %s \n',plotname123)
%%
close all;
s4ContourDiff = s4Contour-s4Contour(1,:);
s4ContourDiff(s4ContourDiff<0) = 0;
pdiffHandle = cell(numel(SysEffectsAll),1-1);
f999 = figure('Units','normalized','Position',[-1,0.5,0.5,0.5]);
for i=1:numel(SysEffectsAll)-1
    if  i==numel(SysEffectsAll)-1
        myColor = rgb('Black');
    else
        myColor = Colors((i+1)*floor(256/numel(SysEffectsAll)),:);
    end
pdiffHandle{i} = plot(s4ContourDiff(i+1,:),m4Contour,...
    'LineWidth',2,'Color',myColor,'LineStyle',LineStyle{i+1});
hold on;
end
set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel(sprintf('|{\\itU}_{e4}|^2(stat. + 1 syst.) - |{\\itU}_{e4}|^2(stat.)'));
ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
PrettyFigureFormat('FontSize',22);
xlim([9.99e-06,0.3])
ylim([min(m4Contour),max(m4Contour)])
leg = legend([pdiffHandle{:}],SysEffectLeg{2:end});
PrettyLegendFormat(leg);
leg.Location = 'east';

%save
plotname999 = sprintf('%sBudget%.0f_SystBreakdownDiff.png',extractBefore(S.DefPlotName,'Budget'),A.SysBudget);
print(f999,plotname999,'-dpng','-r300');
fprintf('save plot to %s \n',plotname999)

% save zoom
ylim([250 400])
xlim([1.5e-05 3.5e-02])
plotname999 = sprintf('%sBudget%.0f_SystBreakdownDiffZoom.png',extractBefore(S.DefPlotName,'Budget'),A.SysBudget);
print(f999,plotname999,'-dpng','-r300');
fprintf('save plot to %s \n',plotname999)
