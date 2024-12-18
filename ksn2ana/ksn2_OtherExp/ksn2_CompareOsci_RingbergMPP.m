% ksn2 calculate chi2 grid search
%% settings that might change
chi2 = 'chi2CMShape';
DataType = 'Twin';


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
Case = 5; % only ksn1
 
if Case==1 || Case==2
    Arg = {'Solar','OFF',...
    'BestGA','OFF',...
    'Best','OFF',...
    'RAA','OFF',...%'RAA_GA','ON',...
    'Stereo','OFF','DANSS','OFF','Prospect','OFF',...
    'DoubleChooz','OFF',...
    'DayaBay','OFF',...%
    'Neutrino4','OFF',...
    'Mainz','OFF','Troitsk','OFF',...
    'NuBetaBeta','OFF',...
    'FinalSensitivity','OFF',...
    'Combi','OFF'};
elseif Case==3
  Arg = {'Solar','OFF',...
    'BestGA','OFF',...
    'Best','OFF',...
    'RAA','OFF',...%'RAA_GA','ON',...
    'Stereo','OFF','DANSS','OFF','Prospect','OFF',...
    'DoubleChooz','OFF',...
    'DayaBay','OFF',...%
    'Neutrino4','OFF',...% 'Mainz','OFF','Troitsk','OFF',...
    'NuBetaBeta','OFF',...
    'FinalSensitivity','ON',...
    'Combi','OFF'};
elseif Case==4
      Arg = {'Solar','OFF',...
    'BestGA','ON',...
    'Best','OFF',...
    'RAA','ON',...%'RAA_GA','ON',...
    'Stereo','OFF','DANSS','OFF','Prospect','OFF',...
    'DoubleChooz','OFF',...
    'DayaBay','OFF',...%
    'Neutrino4','ON',...% 'Mainz','OFF','Troitsk','OFF',...
    'NuBetaBeta','OFF',...
    'FinalSensitivity','ON',...
    'Combi','OFF'};
elseif Case==5
      Arg = {'Solar','OFF',...
    'BestGA','ON',...
    'Best','OFF',...
    'RAA','ON',...%'RAA_GA','ON',...
   ... %'Stereo','OFF','DANSS','OFF','Prospect','OFF',...
   ... %'DoubleChooz','OFF',...
   ... %'DayaBay','OFF',...%
    'Neutrino4','ON',...% 'Mainz','OFF','Troitsk','OFF',... %'NuBetaBeta','OFF',...
    'FinalSensitivity','ON',...
    'Combi','OFF'};
end

[legHandle,legStr,pHandle] = S.ContourPlotOsci(Arg{:},...
     'SavePlot','OFF','Style','PRL','BestFit','ON','LineStyle','-.');

% load combined result for fixed m_nu
savedir_fix = [getenv('SamakPath'),'ksn2ana/ksn2_RunCombination/results/'];
MakeDir(savedir_fix)
savefile_fix = sprintf('%sksn21_Combination_ReAna_Real.mat',savedir_fix);

% load file
if exist(savefile_fix,'file') 
    dfixT = importdata(savefile_fix);
    dfixD = importdata(savefile_fix);
    fprintf('load results from file %s \n',savefile_fix)
else
     fprintf('load results from file %s \n',savefile_fix)
    return
end

mNu4Sq_fix   = dfixT.mNu4Sq_contour_12;
sin2T4Sq_fix  = 4*dfixT.sin2T4_contour_12.*(1-dfixT.sin2T4_contour_12);
mNu4Sq_fix1   = dfixD.mNu4Sq_contour_1;
sin2T4Sq_fix1  = 4*dfixD.sin2T4_contour_1.*(1-dfixD.sin2T4_contour_1);

% include n plot
hold on;
if Case ~=1
    pfix12 = plot(sin2T4Sq_fix,mNu4Sq_fix,'-','Color',rgb('Navy'),'LineWidth',2.5);
end
pfix1 = plot(sin2T4Sq_fix1,mNu4Sq_fix1,':','Color',rgb('SkyBlue'),'LineWidth',2.5);

ax = gca;
set(gca,'FontSize',24);
ax.XLabel.FontSize = 26;
ax.YLabel.FontSize = 26;
  ylabel(sprintf('\\Delta{\\itm}_{41}^2 (eV^{ 2})'));    
            
legStrCombi = legStr;
pNone = plot(NaN,NaN,'LineStyle','none');
% legend
if Case==1
    legHandleAll = [pfix1];
    legHandle{1}.delete;
    legStrCombi = sprintf('KATRIN (KSN1) 95%% C.L.');
   % legStrCombi = {legStrCombi{1:3},sprintf('KATRIN (KSN1) 95%% C.L.'),...
   %     sprintf('KATRIN (KSN1+2) 95%% C.L.')};
elseif Case==2
    legHandleAll = [pfix1,legHandle{1},pfix12];
      legStrCombi ={sprintf('KATRIN (KSN1) 95%% C.L.'),...
          sprintf('KATRIN (KSN2) 95%% C.L.'),...
          sprintf('KATRIN (KSN1+2) 95%% C.L.')};
      ch=get(gca,'children');
      set(gca,'children',[ch(4) ; ch(1:3)]) ;
elseif Case==3
    legHandleAll = [pfix1,legHandle{3},pfix12,legHandle{4},legHandle{5},legHandle{1:2}];
    legStrCombi ={sprintf('KATRIN (KSN1) 95%% C.L.'),...
        sprintf('KATRIN (KSN2) 95%% C.L.'),...
        sprintf('KATRIN (KSN1+2) 95%% C.L.'),...
       legStrCombi{4},legStrCombi{5},...
       legStrCombi{1:2}};
    ch=get(gca,'children');
    set(gca,'children',[ch(6) ; ch(1:5); ch(7:end)]) ;
else 
    legHandleAll = [pfix1,legHandle{end-2},pfix12,legHandle{end-1:end},legHandle{1:end-3}];
    legStrCombi ={sprintf('KATRIN (KSN1) 95%% C.L.'),...
        sprintf('KATRIN (KSN2) 95%% C.L.'),...
        sprintf('KATRIN (KSN1+2) 95%% C.L.'),...
        legStrCombi{end-1:end},...
        legStrCombi{1:end-3}};
    ch=get(gca,'children');
    set(gca,'children',[ch(end-2) ; ch(1:end-3); ch(end-1:end)]) ;  
end

if Case>1
 
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
%
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
ax.XLabel.FontSize = LocalFontSize+2;
ax.YLabel.FontSize = LocalFontSize+2;
xlabel(sprintf('{sin}^2(2\\theta_{14})'));% change x-label
PrettyFigureFormat('FontSize',LocalFontSize);
set(gcf,'Position',[0.1,0.1,0.6,0.5]);
leg.Location = 'northeastoutside';
leg.FontSize = get(gca,'FontSize');
leg.NumColumns = 1;
%leg.Position(2) = 0.31;
ax.Position(4) = 0.75;
ax.Position(3) = 0.4;
ax.Position(2) = 0.22;
%legHandle{1}.Color = rgb('ForestGreen'); % change DANSS color: Grey -> Green

%
%return
pltdir = [getenv('SamakPath'),'ksn2ana/ksn2_OtherExp/plots/'];
MakeDir(pltdir)
% save as pdf
plotname = sprintf('%sOsciContour_Ringberg_%.0f.pdf',pltdir,Case);
export_fig(gcf,plotname);
fprintf('save plot to %s \n',plotname)         
%print(plotname,'-dpng','-r350');

% plotname = sprintf('%s_OsciContour_%.2gCL_CombiBEST_%.0f.png',S.DefPlotName,S.ConfLevel,Case);
% print(gcf,plotname,'-dpng','-r300');
% fprintf('save plot to %s \n',plotname)  

