% ksn2 calculate chi2 grid search
%% settings that might change
chi2 = 'chi2CMShape';
DataType = 'Twin';
Case = 3;

range = 40;
freePar = 'E0 Norm Bkg';
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
  
if Case ==1
    Arg = {'Solar','OFF',...
    'BestGA','OFF',...
    'Best','OFF',...
    'RAA','ON',...%'RAA_GA','ON',...
   'Stereo','OFF','DANSS','OFF','Prospect','OFF',...
   'DoubleChooz','OFF',...
    'DayaBay','OFF',...%
    'Neutrino4','OFF',...% 'Mainz','OFF','Troitsk','OFF',
   'NuBetaBeta','OFF',...
    'FinalSensitivity','OFF',...
    'Combi','OFF'};
elseif Case==2
   Arg = {'Solar','OFF',...
    'BestGA','OFF','Best','OFF',...
    'RAA','ON',...%'RAA_GA','ON',...
    'DoubleChooz','OFF','DayaBay','OFF',...%
    'Neutrino4','OFF',...% 'Mainz','OFF','Troitsk','OFF',
    'NuBetaBeta','OFF',...
    'FinalSensitivity','OFF',...
    'Combi','OFF'};
elseif Case==3
    Arg = {'Solar','OFF',...
        'BestGA','ON','Best','OFF',...
        'RAA','ON','RAA_GA','OFF',...
        'DoubleChooz','OFF','DayaBay','OFF',...%
        'Neutrino4','ON',...% 'Mainz','OFF','Troitsk','OFF',
        'NuBetaBeta','OFF',...
        'FinalSensitivity','ON',...
        'Combi','OFF'};
end

[legHandle,legStr,pHandle] = S.ContourPlotOsci(Arg{:},...
     'SavePlot','OFF','Style','PRL','BestFit','ON');

% load combined result for fixed m_nu
savedir_fix = [getenv('SamakPath'),'ksn2ana/ksn2_RunCombination/results/'];
MakeDir(savedir_fix)
savefile_fixT = sprintf('%sksn21_Combination_ReAna_Twin.mat',savedir_fix);
savefile_fixD = sprintf('%sksn21_Combination_ReAna_Real.mat',savedir_fix);

% load file
if exist(savefile_fixT,'file') 
    dfixT = importdata(savefile_fixT);
    dfixD = importdata(savefile_fixD);
    fprintf('load results from file %s \n',savefile_fixT)
else
     fprintf('load results from file %s \n',savefile_fixT)
    return
end

mNu4Sq_fix   = dfixT.mNu4Sq_contour_12;
sin2T4Sq_fix  = 4*dfixT.sin2T4_contour_12.*(1-dfixT.sin2T4_contour_12);
mNu4Sq_fix1   = dfixD.mNu4Sq_contour_1;
sin2T4Sq_fix1  = 4*dfixD.sin2T4_contour_1.*(1-dfixD.sin2T4_contour_1);

% include n plot
hold on;
pfix12 = plot(sin2T4Sq_fix,mNu4Sq_fix,'-.','Color',rgb('Navy'),'LineWidth',2.5);
pfix1 = plot(sin2T4Sq_fix1,mNu4Sq_fix1,'-','Color',rgb('SkyBlue'),'LineWidth',2.5);

ax = gca;
set(gca,'FontSize',24);
ax.XLabel.FontSize = 26;
ax.YLabel.FontSize = 26;
 
legStrCombi = legStr;
pNone = plot(NaN,NaN,'LineStyle','none');
% legend
if Case==1
    legHandleAll = [legHandle{1:3},pfix1,pfix12,pNone];
    legStrCombi = {legStrCombi{1:3},sprintf('KATRIN (KSN1) 95%% C.L.'),...
        sprintf('KATRIN (KSN1+2)'), 'sensitivity 95% C.L.'};
elseif Case==2
    legHandleAll = [legHandle{1:6},pfix1,pfix12,pNone];
    legStrCombi = {legStrCombi{1:6},sprintf('KATRIN (KSN1) 95%% C.L.'),...
        sprintf('KATRIN (KSN1+2)'), 'sensitivity 95% C.L.'};
elseif Case==3
    legHandleAll = [legHandle{1:8},pfix1,pfix12,pNone,legHandle{end-1:end}];
    legStrCombi = {legStrCombi{1:8},...
        sprintf('KATRIN (KSN1) 95%% C.L.'),...
        sprintf('KATRIN (KSN1+2)'), ...
        'sensitivity 95%C.L.',...
        legStrCombi{end-1:end}};
end


% if strcmp(FinalSensitivity,'ON')
%     legStrCombi{end-2} = sprintf('KATRIN (KSN2) sensitivity 95%% C.L.');
%     
%     legStrCombi = {legStrCombi{1:end-3},...
%         sprintf('KATRIN (KSN1) 95%% C.L.'),...%   [legStr{end-2},' (KSN1)'],...
%         sprintf('KATRIN (KSN2) \nsensitivity 95%% C.L.'),...%   [legStr{end-2},' (KSN2)'],...
%         sprintf('KATRIN (KSN1+2) \nsensitivity 95%% C.L.'),...%    [legStr{end-2},' (KSN1+2)'],...
%         legStrCombi{end-1:end}};
%     legHandleAll = [legHandle{1:end-3},...
%         pfix1,legHandle{end-2},pfix12,...
%         legHandle{end-1},legHandle{end}];
% else
%     legStrCombi{end} = [legStr{end},' (KSN2)'];
%     legStrCombi = {legStrCombi{1:end-1},[legStr{end},' (KSN1)'],legStrCombi{end},[legStr{end},' (KSN1+2)']};
%     legHandleAll = [legHandle{1:end-1},pfix1,legHandle{end},pfix12];
% end

leg = legend(legHandleAll,legStrCombi);
leg.FontSize=19;
leg.Position(1) = 0.005;
leg.Position(2) = 0.7;


ax.Position(4) = 0.55;
xlim([4e-03 1]);
%
LocalFontSize = 20;
set(gca,'FontSize',LocalFontSize);
ax = gca;
ax.XLabel.FontSize = LocalFontSize;
ax.YLabel.FontSize = LocalFontSize;

% modifications for Susanne:
%legHandle{9}.delete; % remove KSN-2 only
if Case==1
legHandle{4}.delete; % remove KSN-2 only
elseif Case==2
    legHandle{7}.delete; 
    legHandle{1}.Color = rgb('ForestGreen');
elseif Case==3
    legHandle{9}.delete;
    legHandle{1}.Color = rgb('ForestGreen');
end

PrettyFigureFormat('FontSize',LocalFontSize);
set(gcf,'Position',[0.1,0.1,0.6,0.5]);
leg.Location = 'northeastoutside';
leg.FontSize = get(gca,'FontSize');
leg.NumColumns = 1;
%leg.Position(2) = 0.31;
ax.Position(4) = 0.75;
ax.Position(3) = 0.4;
ax.Position(2) = 0.22;
xlabel(sprintf('{sin}^2(2\\theta_{14})'));% change x-label
%legHandle{1}.Color = rgb('ForestGreen'); % change DANSS color: Grey -> Green



% save as pdf
plotname = sprintf('%s_OsciContour_%.2gCL_CombiBEST_%.0f.pdf',S.DefPlotName,S.ConfLevel,Case);
export_fig(gcf,plotname);
fprintf('save plot to %s \n',plotname)         
%print('KSN2_Osci.png','-dpng','-r70');

plotname = sprintf('%s_OsciContour_%.2gCL_CombiBEST_%.0f.png',S.DefPlotName,S.ConfLevel,Case);
print(gcf,plotname,'-dpng','-r300');
fprintf('save plot to %s \n',plotname)  

