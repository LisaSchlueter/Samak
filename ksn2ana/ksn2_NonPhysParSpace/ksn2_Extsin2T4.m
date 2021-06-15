% ksn2
% m2 nuisance parameter unconstrained
% closer look at area [0.5,1] and compare
%% settings that might change
chi2 = 'chi2Stat';
DataType = 'Twin';
nGridSteps = 30;
range = 40;
Mode = 'Display';
%% configure RunAnalysis object
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
    'range',range};

S = SterileAnalysis(SterileArg{:});

if strcmp(A.DataSet,'Twin')
    ExtmNu4Sq = 'OFF';
else
    ExtmNu4Sq = 'ON';
end
%%
switch Mode
    case 'Compute'
        S.GridSearch('ExtmNu4Sq',ExtmNu4Sq,'mNu4SqTestGrid',5); 
        A = MultiRunAnalysis(RunAnaArg{:});
        S = SterileAnalysis(SterileArg{:});    
        S.GridSearch('ExtmNu4Sq',ExtmNu4Sq,'mNu4SqTestGrid',5.5);
    case 'Display'
        S.LoadGridFile('ExtmNu4Sq',ExtmNu4Sq,'mNu4SqTestGrid',5,'CheckLargerN','ON');
        S.InterpMode = 'spline';
         S.InterpMode = 'lin';
        mNu4Sq1 =  S.mNu4Sq;
        sin2T41 =  S.sin2T4;
        chi21 =  S.chi2;
        mNuSq1 =  S.mNuSq;
        E01   = S.E0;
        S.Interp1Grid;
        PlotPar1 = S.mNuSq;
        
        % isolines
        IsoVec = [-10,-5,-1,0,0.3,1.1,2,10,20];
        GetFigure;
        
        for i=1:numel(ContourVec)
            [~,p1] = contour3(S.sin2T4,S.mNu4Sq,PlotPar1,[IsoVec(i),IsoVec(i)],...
                'Color',rgb('Silver'),'ShowText','on','LineWidth',1.5,'LabelSpacing',380);
            hold on;
        end
        view(2)
        grid off
        preg = S.ContourPlot('HoldON','ON');
        
     %  
        S.LoadGridFile('ExtmNu4Sq',ExtmNu4Sq,'mNu4SqTestGrid',5.5);
        mNu4Sq2 =  S.mNu4Sq;
        mNuSq2 =  S.mNuSq;
        E02   = S.E0;
        sin2T42 =  S.sin2T4;
        chi22 =  S.chi2;
        %
        S.Interp1Grid;%('Maxm4Sq',10^2);
          PlotPar2 = S.mNuSq;
        for i=1:numel(ContourVec)
            [~,p1] = contour3(S.sin2T4,S.mNu4Sq,PlotPar2,[IsoVec(i),IsoVec(i)],...
                'Color',rgb('Silver'),'ShowText','on','LineWidth',1.5,'LabelSpacing',380);
            hold on;
        end
        pext= S.ContourPlot('HoldOn','ON','LineStyle','-.','Color',rgb('Orange'));
        set(gca,'XScale','lin')
        xlim([0 1])
        
        return
        %% merge
        S.InterpMode = 'lin';
        S.mNu4Sq = repmat(mNu4Sq1(:,1),1,S.nGridSteps*2-1);
        S.sin2T4 = repmat([sin2T41(1,:),sin2T42(1,2:end)],S.nGridSteps,1);
        S.chi2  = [chi21',chi22(2:end,:)']';
        S.mNuSq = [mNuSq1',mNuSq2(2:end,:)']';
     
        S.E0 = [E01',E02(2:end,:)']';
        % chi2g(chi2g>5.99)=NaN;
        %  S.GridPlot;
        S.Interp1Grid
            S.mNuSq(abs(S.mNuSq)>10) = NaN;
            S.GridPlotFitPar
end


