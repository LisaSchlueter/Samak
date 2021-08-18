% ksn2: all systematics minus 1 effects
MinusSysEffect = 'allmLongPlasmaPT';%'BkgPT';

savedir = [getenv('SamakPath'),'ksn2ana/ksn2_Systematics/results/'];
savename = sprintf('%sksn2_StatOverSyst_Sysm1_%s.mat',savedir,MinusSysEffect);
%savename = sprintf('%sksn2_StatOverSyst_Top3.mat',savedir);

if exist(savename,'file')
else
    %% settings that might change
    chi2 = 'chi2Stat';
    DataType = 'Twin';
    nGridSteps = 30;
    range = 40;
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
        'LoadGridArg',LoadGridArg}
    
%     %%
     S = SterileAnalysis(SterileArg{:});
    S.RunAnaObj.chi2 = 'chi2Stat';
    S.RunAnaObj.NonPoissonScaleFactor = 1;
    S.LoadGridFile(S.LoadGridArg{:});
    S.Interp1Grid('MaxM4Sq',38^2);
    S.ContourPlot;
    sin2T4_contour_stat = S.sin2T4_contour;
    mNu4Sq_contour_stat = S.mNu4Sq_contour;
    
    %%
    S.RunAnaObj.chi2 = 'chi2CMShape';
    if strcmp(MinusSysEffect,'NP')
        S.RunAnaObj.NonPoissonScaleFactor = 1;
        S.SysEffect = 'all';
    elseif strcmp(MinusSysEffect,'BkgPT')
        S.RunAnaObj.NonPoissonScaleFactor = 1.112;
        S.SysEffect = 'allmBkgPT';
    elseif strcmp(MinusSysEffect,'allmLongPlasma')
         S.RunAnaObj.NonPoissonScaleFactor = 1.112;
        S.SysEffect = 'allmLongPlasma';  
    elseif  strcmp(MinusSysEffect,'allmLongPlasmaPT')
          S.SysEffect = 'allmLongPlasmaPT';  
        strcmp(MinusSysEffect,'allmLongPlasmaPT');
    end
  %  S.SysEffect = 'KSN2Top3';
    
%      S.GridSearch(S.LoadGridArg{:});
%      return
    S.LoadGridFile(S.LoadGridArg{:});
     S.Interp1Grid('MaxM4Sq',38^2);
    S.ContourPlot('HoldOn','ON');
    sin2T4_contour_cm = S.sin2T4_contour;
    mNu4Sq_contour_cm = S.mNu4Sq_contour;
    
    %%
     mNu4Sq    = linspace(max([min(mNu4Sq_contour_stat),min(mNu4Sq_contour_cm)]),min([max(mNu4Sq_contour_stat),max(mNu4Sq_contour_cm)]),1e3);
    sin2T4_cm = interp1(mNu4Sq_contour_cm,sin2T4_contour_cm,mNu4Sq,'spline');
    sin2T4_stat = interp1(mNu4Sq_contour_stat,sin2T4_contour_stat,mNu4Sq,'spline');
    sin2T4_sys = sqrt(sin2T4_cm.^2-sin2T4_stat.^2);
    
    save(savename,'mNu4Sq','sin2T4_cm','sin2T4_stat','sin2T4_sys')
   fprintf('save to %s \n',savename);
    
end