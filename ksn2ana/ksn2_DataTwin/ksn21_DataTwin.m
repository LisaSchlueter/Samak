% plot exclusion and sensitivity for ksn1 and ksn2 (and Combination)
Combi = 'OFF';
chi2 = 'chi2CMShape';
freePar = 'mNu E0 Norm Bkg';

savedir = [getenv('SamakPath'),'ksn2ana/ksn2_DataTwin/results/'];
MakeDir(savedir);
savefile2 = sprintf('%sksn2_DataTwinContour_%s_%s.mat',savedir,chi2,strrep(freePar,' ',''));
if exist(savefile2,'file') 
    load(savefile2);
    fprintf('load from file %s \n',savefile2);
else  
    DataType = 'Real';
    nGridSteps = 40;
    range = 40;
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
        'nGridSteps',nGridSteps,...
        'SmartGrid','OFF',...
        'RecomputeFlag','OFF',...
        'SysEffect','all',...
        'RandMC','OFF',...
        'range',range,...
        'LoadGridArg',{'mNu4SqTestGrid',5}};
    
    S = SterileAnalysis(SterileArg{:});
    
    % twins
    S.RunAnaObj.DataType = 'Twin';
    S.LoadGridFile('CheckSmallerN','ON',S.LoadGridArg{:},'ExtmNu4Sq','OFF');
    S.Interp1Grid('RecomputeFlag','ON','MinM4Sq',1,'MaxM4Sq',40^2);
    S.ContourPlot('BestFit','ON'); close;
    if S.chi2_Null<S.chi2_ref
        chi2_T2 = S.chi2-S.chi2_Null;
    else
        chi2_T2 = S.chi2-S.chi2_ref;
    end
    sin2T4_contourT2 = S.sin2T4_contour;
    mNu4Sq_contourT2 = S.mNu4Sq_contour;
    % data
    S.RunAnaObj.DataType = 'Real';
    if ~contains(freePar,'mNu')
        S.LoadGridFile('CheckSmallerN','ON',S.LoadGridArg{:},'ExtmNu4Sq','ON','Extsin2T4','ON');
        S.InterpMode = 'Mix';
        S.Interp1Grid('RecomputeFlag','ON');
    else
        S.LoadGridFile('CheckSmallerN','ON',S.LoadGridArg{:},'ExtmNu4Sq','ON');
        S.Interp1Grid('RecomputeFlag','ON','MinM4Sq',1,'MaxM4Sq',40^2);
    end
    
    S.ContourPlot('BestFit','ON'); close;
    sin2T4_contourD2 = S.sin2T4_contour;
    mNu4Sq_contourD2 = S.mNu4Sq_contour;
    sin2T4_bfD2 = S.sin2T4_bf;
    mNu4Sq_bfD2 = S.mNu4Sq_bf;
    chi2min_bfD2 = S.chi2_bf;
    
    if ~contains(freePar,'mNu')
        S.LoadGridFile('CheckSmallerN','ON',S.LoadGridArg{:},'ExtmNu4Sq','ON','Extsin2T4','OFF');
        S.InterpMode = 'spline';
        S.Interp1Grid('RecomputeFlag','ON','MinM4Sq',1,'MaxM4Sq',40^2);
        chi2_D2 = S.chi2-chi2min_bfD2;
    else
        chi2_D2 = S.chi2-chi2min_bfD2;
    end
    
    save(savefile2,'chi2_T2','chi2_D2','sin2T4_contourT2','sin2T4_contourD2',...
        'mNu4Sq_contourT2','mNu4Sq_contourD2','mNu4Sq_bfD2','sin2T4_bfD2','chi2min_bfD2',...
        'RunAnaArg','SterileArg')
end


savefile1 = sprintf('%sksn1_DataTwinContour_%s_%s.mat',savedir,chi2,strrep(freePar,' ',''));
if exist(savefile1,'file')
    load(savefile1);
      fprintf('load from file %s \n',savefile1);
else
    %% ksn1 (re-ana)
    nGridSteps            = 50;
    DataType              = 'Twin';
    % configure RunAnalysis object
    if strcmp(chi2,'chi2Stat')
        NonPoissonScaleFactor = 1;
    elseif  strcmp(chi2,'chi2CMShape')
        NonPoissonScaleFactor = 1.064;
    end
    range = 40;
    Real = MultiRunAnalysis('RunList','KNM1',...
        'chi2',chi2,...
        'DataType',DataType,...
        'fixPar',freePar,...
        'NonPoissonScaleFactor',NonPoissonScaleFactor,...
        'SysBudget',200,...
        'minuitOpt','min ; minos',...
        'FSDFlag','KNM2_0p1eV',...
        'ELossFlag','KatrinT2A20',...
        'AngularTFFlag','ON',...
        'SynchrotronFlag','ON',...
        'RadiativeFlag','ON',...
        'DopplerEffectFlag','FSD',...
        'BKG_PtSlope',-2.2*1e-06);
    Real.exclDataStart = Real.GetexclDataStart(40);
    % configure Sterile analysis object
    SterileArg = {'RunAnaObj',Real,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
        'nGridSteps',nGridSteps,...
        'SmartGrid','OFF',...
        'RecomputeFlag','OFF',...
        'SysEffect','all',...
        'RandMC','OFF',...
        'range',range};
    %%
    S = SterileAnalysis(SterileArg{:});

    S.LoadGridFile('ExtmNu4Sq','OFF','mNu4SqTestGrid',2);
    S.Interp1Grid('RecomputeFlag','ON','Maxm4Sq',37^2);
    S.ContourPlot('BestFit','ON'); close;

     if  contains(freePar,'mNu')
    sin2T4_contourT1 = S.sin2T4_contour(:,3)';
    mNu4Sq_contourT1 = S.mNu4Sq_contour(:,3)';
     else
            sin2T4_contourT1 = S.sin2T4_contour;
    mNu4Sq_contourT1 = S.mNu4Sq_contour;
     end
    %%
    S.RunAnaObj.DataType = 'Real';
    
    if  contains(freePar,'mNu')
        S.InterpMode = 'Mix';
    end
    S.LoadGridFile('ExtmNu4Sq','OFF','mNu4SqTestGrid',2);
    S.Interp1Grid;
    S.ContourPlot('BestFit','ON'); close;
    sin2T4_contourD1 = S.sin2T4_contour;
    mNu4Sq_contourD1 = S.mNu4Sq_contour;
    sin2T4_bfD1 = S.sin2T4_bf;
    mNu4Sq_bfD1 = S.mNu4Sq_bf;
    chi2min_bfD1 = S.chi2_bf;
    
    save(savefile1,'sin2T4_contourT1','sin2T4_contourD1',...
                  'mNu4Sq_contourT1','mNu4Sq_contourD1',...
                 'mNu4Sq_bfD1','sin2T4_bfD1','chi2min_bfD1')
end

%% load combinatinom
if strcmp(Combi,'ON') && ~contains(freePar,'mNu')
    fileCT = sprintf('%sksn2ana/ksn2_RunCombination/results/ksn21_Combination_ReAna_Twin.mat',getenv('SamakPath'));
    dCT = importdata(fileCT);
    fileCD = sprintf('%sksn2ana/ksn2_RunCombination/results/ksn21_Combination_ReAna_Real.mat',getenv('SamakPath'));
    dCD = importdata(fileCD);
    fprintf('load from file %s \n',fileCT);
    fprintf('load from file %s \n',fileCD);
    mNu4Sq_CT    = dCT.mNu4Sq_contour_12;
    sin2T4_CT    = dCT.sin2T4_contour_12;
    mNu4Sq_CD    = dCD.mNu4Sq_contour_12;
    sin2T4_CD    = dCD.sin2T4_contour_12;
    mNu4Sq_CD_bf = dCD.mNu4Sq_bf_12;
    sin2T4_CD_bf = dCD.sin2T4_bf_12;
elseif strcmp(Combi,'ON') && contains(freePar,'mNu')
    CTdir = [getenv('SamakPath'),'SterileAnalysis/GridSearchFiles/Combi/'];
    fileCT = sprintf('%sTwin/KSN12Combi_ReAna_GridSearch_Twin_mNuE0BkgNorm_Uniform_chi2CMShape_30nGrid.mat',CTdir);
    dCT = importdata(fileCT);
    fileCD = sprintf('%sReal/KSN12Combi_ReAna_GridSearch_Real_mNuE0BkgNorm_Uniform_chi2CMShape_30nGrid.mat',CTdir);
    dCD = importdata(fileCD);
    fprintf('load from file %s \n',fileCT);
    fprintf('load from file %s \n',fileCD);
    mNu4Sq_CT    = dCT.mNu4Sq_contour;
    sin2T4_CT    = dCT.sin2T4_contour;
    mNu4Sq_CD    = dCD.mNu4Sq_contour;
    sin2T4_CD    = dCD.sin2T4_contour;
    mNu4Sq_CD_bf = dCD.mNu4Sq_bf;
    sin2T4_CD_bf = dCD.sin2T4_bf;
end
%% plot
GetFigure;
pD1 = plot(sin2T4_contourD1,mNu4Sq_contourD1,'-','LineWidth',2,'Color',rgb('DodgerBlue'));
hold on;
pT1 = plot(sin2T4_contourT1,mNu4Sq_contourT1,':','LineWidth',2,'Color',pD1.Color);
pD2 = plot(sin2T4_contourD2,mNu4Sq_contourD2,'-','LineWidth',2,'Color',rgb('Orange'));
pT2 = plot(sin2T4_contourT2,mNu4Sq_contourT2,':','LineWidth',2,'Color',pD2.Color);
pD1bf = plot(sin2T4_bfD1,mNu4Sq_bfD1,'x','Color',pD1.Color,'LineWidth',2,'MarkerSize',9);
pD2bf = plot(sin2T4_bfD2,mNu4Sq_bfD2,'o','Color',pD2.Color,'LineWidth',2,'MarkerSize',pD1bf.MarkerSize);
set(gca,'xScale','log');
set(gca,'yScale','log');
PrettyFigureFormat('FontSize',22)
ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
xlabel(sprintf('|{\\itU}_{e4}|^2'));
xlim([3e-03 0.5]);
ylim([1 1800]);
if strcmp(Combi,'ON') 
    pCT = plot(sin2T4_CT,mNu4Sq_CT,':','LineWidth',2,'Color',rgb('Black'));
    pCD = plot(sin2T4_CD,mNu4Sq_CD,'-','LineWidth',2,'Color',pCT.Color);
    pCDbf = plot(sin2T4_CD_bf,mNu4Sq_CD_bf,'*','MarkerSize',pD1bf.MarkerSize,'LineWidth',2,'Color',pCD.Color);
    leg = legend([pT1,pT2,pCT,pD1,pD2,pCD],...
        'KSN-1 twin', 'KSN-2 twin','KSN-1+2 twin',...
        'KSN-1 data','KSN-2 data','KSN-1+2 data',...
        'Location','southwest');
else
    leg = legend([pT1,pT2,pD1,pD2],'KSN-1 twin','KSN-2 twin','KSN-1 data','KSN-2 data','Location','southwest');
end
leg.NumColumns = 2;
PrettyLegendFormat(leg);
if  contains(freePar,'mNu')
    leg.Title.String =(sprintf('{\\itm}_\\nu free')); 
else
    leg.Title.String =(sprintf('{\\itm}_\\nu fix'));   
end
leg.Title.FontWeight = 'normal';
leg.Title.FontSize = get(gca,'FontSize');

pltdir = strrep(savedir,'results','plots');
MakeDir(pltdir);
pltname =  sprintf('%sksn21_DataTwinContour_%s_%s_Combi%s.png',pltdir,chi2,strrep(freePar,' ',''),Combi);
%print(gcf,pltname,'-dpng','-r350');
fprintf('save plot to %s \n',pltname);