% ksn2
% m2 nuisance parameter unconstrained
% closer look at area [0.5,1] and compare
% no interpolation
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
        
        % grid 1
        S.LoadGridFile('ExtmNu4Sq',ExtmNu4Sq,'mNu4SqTestGrid',5,'CheckLargerN','ON');
        mNu4Sq1  =  S.mNu4Sq;
        sin2T41  =  S.sin2T4;
        chi21    =  S.chi2;
        mNuSq1   =  S.mNuSq;
        E01      = S.E0;
        PlotPar1 = S.mNuSq;
        
        % grid 2
        S.LoadGridFile('ExtmNu4Sq',ExtmNu4Sq,'mNu4SqTestGrid',5.5);
        mNu4Sq2 =  S.mNu4Sq;
        mNuSq2  =  S.mNuSq;
        E02     = S.E0;
        sin2T42 =  S.sin2T4;
        chi22   =  S.chi2;
        
        %% merge
        S.mNu4Sq = repmat(mNu4Sq1(:,1),1,S.nGridSteps*2-1);
        S.sin2T4 = repmat([sin2T41(1,:),sin2T42(1,2:end)],S.nGridSteps,1);
        S.chi2  = [chi21',chi22(2:end,:)']';
        S.mNuSq = [mNuSq1',mNuSq2(2:end,:)']';
        S.E0 = [E01',E02(2:end,:)']';
        
        S.InterpMode = 'spline';
        S.Interp1Grid('Maxm4Sq',38^2)
        
        
        
        f1 = figure('Units','normalized','Position',[0.1,0.1,0.7,0.55]);
        % subplot 2
        subplot(1,2,2)
        
        
        ContourVec = [-700 -300 -100 -50 -10 -3];
        S.IsoPlotFitPar('HoldOn','ON','ContourVec',ContourVec);
        title('')
        legKeep = get(gca,'Legend');
        
        set(gca,'XScale','lin')
        xlim([0.5 1])
        ylabel(''); xlabel('');
        ax1 = gca;
        yticklabels('');
        
        
        % only first grid
        subplot(1,2,1)
        
        S.mNu4Sq = mNu4Sq1;
        S.sin2T4 = sin2T41;
        S.chi2  = chi21;
        S.mNuSq = mNuSq1;
        S.E0 = E01;
        S.InterpMode = 'spline';
        S.Interp1Grid('Maxm4Sq',38^2)
        
        S.IsoPlotFitPar('HoldOn','ON','ContourVec',[-10 -3 -1 -0.1 0 0.1 0.3 1 5 10]);
        leg = get(gca,'Legend');
        leg.delete;
        title('')
        
        ax2 = gca;
        ax2.Position(3) = 0.35;  ax1.Position(3) =  ax2.Position(3);
        ax2.Position(4) = 0.65;  ax1.Position(4) =  ax2.Position(4);
        ax2.Position(1) = 0.1;   ax1.Position(1) = 0.452;
        ax2.Position(2) = 0.16;  ax1.Position(2) = ax2.Position(2);
        ax2.XLabel.Position(1) = 0.5;
        legKeep.Position(1) = 0.36;
        legKeep.Position(2) = 0.83;
        %%
        pltname = [S.DefPlotName,sprintf('Extsin2T4_GridPlt_Txt%s_NegmNuSq%s.png',ContourTxt,NegmNuSq)];
        print(pltname,'-dpng','-r350');
end

        


