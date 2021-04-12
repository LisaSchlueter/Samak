% ksn2 calculate chi2 grid search
% compare m2 free, m2 nuisance parameter
%% settings that might change
chi2 = 'chi2Stat';
DataType = 'Twin';
nGridSteps = 25;
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
    'FSDFlag','KNM2_0p1eV',...
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
    'range',range};

%%
S = SterileAnalysis(SterileArg{:});
S.InterpMode = 'lin';
S.RunAnaObj.chi2 = 'chi2CMShape';
S.RunAnaObj.NonPoissonScaleFactor = 1;
SysEffectsAll   = {'Stat','FSD','BkgPT','Bkg','NP','LongPlasma','TASR','RF_EL',...
    'RF_BF','RF_RX','FPDeff','Stack','TCoff_OTHER'}; 

HoldOn = 'OFF';
pHandle = cell(numel(SysEffectsAll),1);
LineStyle = {'-','-.',':','--','-','-.',':','--','-','-.',':','--','-.',':','--','-.',':','--','-.',':','--'};
Colors = colormap('jet');

Contour_x = cell(numel(SysEffectsAll),1);
Contour_y = cell(numel(SysEffectsAll),1);

m4Contour = logspace(log10(1.5),log10(34^2),1e3);
s4Contour = zeros(numel(SysEffectsAll),1e3);
for i=1:numel(SysEffectsAll)
    S.SysEffect = SysEffectsAll{i};
    
    if strcmp(SysEffectsAll{i},'NP') || strcmp(SysEffectsAll{i},'Stat')
        S.RunAnaObj.chi2 = 'chi2Stat';
        if strcmp(SysEffectsAll{i},'NP')
            S.RunAnaObj.NonPoissonScaleFactor = 1.112;
        end
      %  S.GridSearch;
        
        S.LoadGridFile;
        S.Interp1Grid('Maxm4Sq',35^2);
        
        S.RunAnaObj.chi2 = 'chi2CMShape';
        S.RunAnaObj.NonPoissonScaleFactor = 1;

    else
     %   S.GridSearch;
        
        S.LoadGridFile;
        S.Interp1Grid('Maxm4Sq',35^2);
    end

    pHandle{i} =  S.ContourPlot('HoldOn',HoldOn,'SavePlot','OFF','Color',Colors(i*floor(256/numel(SysEffectsAll)),:),'LineStyle',LineStyle{i});
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
SysEffectLeg    = {'Stat. only';'Final-state distribution';'Background PT';'Background slope'; 'Non-Poisson background';...
                    'Long. source potential'; 'Tritium activity fluctuations';'Energy-loss function';...
                     'Magnetic fields';'Source scattering';'Detector efficiency';'HV fluctuations';'Theoretical corrections';};
leg = legend([pHandle{:}],SysEffectLeg);
PrettyLegendFormat(leg);
leg.Title.String = 'Stat. + 1 syst.';
leg.Title.FontWeight ='normal';

%save
plotname123 = sprintf('%sBudget%.0f_SystBreakdown.png',extractBefore(S.DefPlotName,'Budget'),A.SysBudget);
print(gcf,plotname123,'-dpng','-r300');
fprintf('save plot to %s \n',plotname123)
%%
close all;
s4ContourDiff = s4Contour-s4Contour(1,:);
s4ContourDiff(s4ContourDiff<0) = 0;
pdiffHandle = cell(numel(SysEffectsAll),1-1);
f999 = figure('Units','normalized','Position',[-1,0.5,0.5,0.5]);
for i=1:numel(SysEffectsAll)-1
pdiffHandle{i} = plot(s4ContourDiff(i+1,:),m4Contour,...
    'LineWidth',2,'Color',Colors((i+1)*floor(256/numel(SysEffectsAll)),:),'LineStyle',LineStyle{i+1});
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

%save
plotname999 = sprintf('%sBudget%.0f_SystBreakdownDiff.png',extractBefore(S.DefPlotName,'Budget'),A.SysBudget);
print(f999,plotname999,'-dpng','-r300');
fprintf('save plot to %s \n',plotname999)