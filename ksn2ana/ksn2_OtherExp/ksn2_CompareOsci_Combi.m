% ksn2 calculate chi2 grid search
%% settings that might change
chi2 = 'chi2CMShape';
DataType = 'Real';

range = 40;
freePar = 'E0 Norm Bkg';
LegLog = 2;

%% configure RunAnalysis object
if strcmp(chi2,'chi2Stat')
    NonPoissonScaleFactor = 1;
elseif  strcmp(chi2,'chi2CMShape')
    NonPoissonScaleFactor = 1.112;
end
RunAnaArg = {'RunList','KNM2_Prompt',...
    'DataType',DataType,...
    'fixPar',freePar,...%free par
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
    'SmartGrid','OFF',...
    'RecomputeFlag','OFF',...
    'SysEffect','all',...
    'RandMC','OFF',...
    'range',range};

%%
S = SterileAnalysis(SterileArg{:});
if strcmp(DataType,'Real')
    S.nGridSteps = 40;
    S.LoadGridArg = {'ExtmNu4Sq','ON','mNu4SqTestGrid',5};
else
     S.LoadGridArg = {'ExtmNu4Sq','ON','mNu4SqTestGrid',2};
     S.nGridSteps = 30;
end
S.LoadGridFile(S.LoadGridArg{:});
S.Interp1Grid;

FinalSensitivity = 'ON';
[legHandle,legStr,pHandle] = S.ContourPlotOsci(....
    'Solar','OFF',...
    'BestGA','ON',...%'Best','ON',...
    'RAA','ON',...%'RAA_GA','ON',...
    'SavePlot','OFF','Style','PRL','BestFit','ON',...
   ...% 'Stereo','OFF','DANSS','OFF','Prospect','OFF','DoubleChooz','OFF',...
   ... % 'DayaBay','OFF','Neutrino4','OFF', 'Mainz','OFF','Troitsk','OFF', 'NuBetaBeta','OFF',...
    'FinalSensitivity',FinalSensitivity,...
    'Combi','OFF',...
    'LineStyle','-.');

% load combined result for fixed m_nu
savedir_fix = [getenv('SamakPath'),'ksn2ana/ksn2_RunCombination/results/'];
MakeDir(savedir_fix)
savefile_fix = sprintf('%sksn21_Combination_ReAna_%s.mat',savedir_fix,DataType);

% load file
if exist(savefile_fix,'file') 
    dfix = importdata(savefile_fix);
    fprintf('load results from file %s \n',savefile_fix)
else
     fprintf('load results from file %s \n',savefile_fix)
    return
end
mNu4Sq_fix   = dfix.mNu4Sq_contour_12;
sin2T4Sq_fix  = 4*dfix.sin2T4_contour_12.*(1-dfix.sin2T4_contour_12);
mNu4Sq_fix1   = dfix.mNu4Sq_contour_1;
sin2T4Sq_fix1  = 4*dfix.sin2T4_contour_1.*(1-dfix.sin2T4_contour_1);

% include n plot
hold on;
pfix12 = plot(sin2T4Sq_fix,mNu4Sq_fix,'-','Color',rgb('Navy'),'LineWidth',3);
pfix1 = plot(sin2T4Sq_fix1,mNu4Sq_fix1,':','Color',rgb('SkyBlue'),'LineWidth',3);

ax = gca;
set(gca,'FontSize',24);
ax.XLabel.FontSize = 26;
ax.YLabel.FontSize = 26;

%% legend
mNuStr = sprintf('{\\itm}_\\nu^2 = 0 eV^2');
legStrCombi = legStr;
if strcmp(FinalSensitivity,'ON')
    legStrCombi{end-1} = [legStr{end-1},' sensitivity'];
    legStrCombi{end} = strrep(legStr{end},'sensitivity',sprintf('(%s)',mNuStr));
    legStrCombi = {legStrCombi{1:end-3},...
       sprintf('KATRIN (KNM1, %s) 95%% C.L.',mNuStr),...% [legStr{end-2},' (KSN1)'],...
       sprintf('KATRIN (KNM2, %s) 95%% C.L.',mNuStr),...% [legStr{end-2},' (KSN2)'],...
      sprintf('KATRIN (KNM1+2, %s) 95%% C.L.',mNuStr),...%  [legStr{end-2},' (KSN1+2)'],...
        legStrCombi{end-1:end}};
    legHandleAll = [legHandle{1:end-3},...
        pfix1,legHandle{end-2},pfix12,...
        legHandle{end-1},legHandle{end}];
else
    legStrCombi{end} = [legStr{end},' (KNM2)'];
    legStrCombi = {legStrCombi{1:end-1},[legStr{end},' (KNM1)'],legStrCombi{end},[legStr{end},' (KNM1+2)']};
    legHandleAll = [legHandle{1:end-1},pfix1,legHandle{end},pfix12];
end

leg = legend(legHandleAll,legStrCombi);
leg.FontSize=19;
leg.Position(1) = 0.005;
leg.Position(2) = 0.7;
ax.Position(4) = 0.55;
xlim([4e-03 1]);

% but ksn2-only in front of ksn1+2
ch=get(gca,'children');
set(gca,'children',[ch(5) ;ch(1:4) ;ch(6:end)]) ;

LocalFontSize = 20;
ax = gca;

% different legend position and style
if LegLog==2
   leg.Location = 'eastoutside'; 
   f1 = gcf;
   f1.Position = [0 0.1 0.6 0.6];
   leg.NumColumns = 1;
   ax.Position(2) = 0.2;
   ax.Position(4) = 0.65;
   ax.Position(3) = 0.45;
   LocalFontSize = 12;
   leg.FontSize = LocalFontSize;
end

set(gca,'FontSize',LocalFontSize);
ax.XLabel.FontSize = LocalFontSize;
ax.YLabel.FontSize = LocalFontSize;


%% save as pdf
plotname = sprintf('%s_OsciContour_%.2gCL_CombiBEST.pdf',S.DefPlotName,S.ConfLevel);
export_fig(gcf,plotname);
fprintf('save plot to %s \n',plotname)         
print('KSN2_Osci.png','-dpng','-r100');



