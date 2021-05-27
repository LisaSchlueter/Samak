% ksn2 grid plot with free nu-mass and isolines
%% settings that might change
chi2 = 'chi2CMShape';
DataType = 'Real';
nGridSteps = 40;
range = 40;
pullFlag = 99;

GridPlot = 'ON';
ContourPlot = 'ON';
%% configure RunAnalysis object
if strcmp(chi2,'chi2Stat')
    NonPoissonScaleFactor = 1;
elseif  strcmp(chi2,'chi2CMShape')
    NonPoissonScaleFactor = 1.112;
end
RunAnaArg = {'RunList','KNM2_Prompt',...
    'DataType',DataType,...
    'fixPar','mNu E0 Norm Bkg',...%free par
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
    'PullFlag',pullFlag,...;%99 = no pull
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
    'LoadGridArg',{'ExtmNu4Sq','ON','mNu4SqTestGrid',5}};

%%
S = SterileAnalysis(SterileArg{:});
S.LoadGridFile(S.LoadGridArg{:});

if strcmp(A.DataType,'Real')
    S.InterpMode = 'spline';
    BF = 'ON';
    % Find Best Fit
    S.LoadGridFile(S.LoadGridArg{:});
    S.Interp1Grid('RecomputeFlag','ON','Maxm4Sq',40^2);
    S.FindBestFit;
    S.FindBestFit('Mode','Imp');
    fprintf('Imp. Best fit: sin2T4 = %.3f m4Sq = %.2feV^2 mNuSq = %.2feV^2 chi2min = %.2f (%.0f dof) \n',...
        S.sin2T4_bf,S.mNu4Sq_bf,S.mNuSq_bf,S.chi2_bf,S.dof);
else
    S.InterpMode = 'spline';
    BF = 'OFF';
end

S.LoadGridFile(S.LoadGridArg{:});
S.Interp1Grid('RecomputeFlag','ON','Maxm4Sq',40^2);
if strcmp(GridPlot,'ON')
    S.GridPlot('BestFit',BF,'Contour','ON','SavePlot','png')
end
%%
if strcmp(ContourPlot,'ON')
    PlotPar = S.mNuSq;
    
    GetFigure;
    
    [~,p1] = contour3(S.sin2T4,S.mNu4Sq,PlotPar,[-1,-1],'Color',rgb('Silver'),'ShowText','on','LineWidth',2);
    hold on;
 %  [~,p1] = contour3(S.sin2T4,S.mNu4Sq,PlotPar,[-0.5,-0.5],'Color',rgb('LimeGreen'),'ShowText','on','LineWidth',1.5);  
    [~,~] = contour3(S.sin2T4,S.mNu4Sq,PlotPar,[1,1],'Color',rgb('Silver'),'ShowText','on','LineWidth',p1.LineWidth);
    [~,~] = contour3(S.sin2T4,S.mNu4Sq,PlotPar,[2,2],'Color',rgb('Silver'),'ShowText','on','LineWidth',p1.LineWidth);
    %[~,~] = contour3(S.sin2T4,S.mNu4Sq,PlotPar,[5,5],'Color',rgb('Silver'),'ShowText','on','LineWidth',p1.LineWidth);
    [~,~] = contour3(S.sin2T4,S.mNu4Sq,PlotPar,[10,10],'Color',rgb('Silver'),'ShowText','on','LineWidth',p1.LineWidth);
    [~,~] = contour3(S.sin2T4,S.mNu4Sq,PlotPar,[0,0],'Color',rgb('Silver'),'ShowText','on','LineWidth',p1.LineWidth);
    [~,~] = contour3(S.sin2T4,S.mNu4Sq,PlotPar,[0.3,0.3],'Color',rgb('Silver'),'ShowText','on','LineWidth',p1.LineWidth);
    
    pcontour = S.ContourPlot('HoldOn','ON','BestFit',BF,'SavePlot','OFF','CL',95,'Color',rgb('ForestGreen'));
    pcontour.LineWidth = 3;
    view(2);
    grid off;
    ylim([0.1 40^2]);
    xlim([1e-03 0.5]);
    % save
    leg = legend([p1,pcontour],...
        sprintf('{\\itm}_\\nu^2 best fit isoline'),...
        sprintf('Free {\\itm}_\\nu^2 unconstrained'),...
        'EdgeColor',rgb('Silver'),'Location','southwest');
    PrettyLegendFormat(leg);
    % save
    name_i = strrep(S.DefPlotName,'_mNuE0BkgNorm','');
    plotname = sprintf('%s_mNuSqNuisance_%.2gCL.png',name_i,S.ConfLevel);
    print(gcf,plotname,'-dpng','-r450');
    fprintf('save plot to %s \n',plotname);
end

%%

