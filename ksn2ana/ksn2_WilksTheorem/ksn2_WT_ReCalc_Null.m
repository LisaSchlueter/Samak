% recalculate chi2 for MC truth
% something weird must have happened...

% grid search on randomized twins
% ksn2 calculae chi2 grid search

Hypothesis = 'H0';
switch Hypothesis
    case 'H0'
        randMC = 1:1e3;
        Twin_sin2T4 = 0;
        Twin_mNu4Sq = 0;
        chi2 = 'chi2CMShape';
    case 'H1'
        randMC = [1:340,384:794];
        Twin_sin2T4 = 0.0240;
        Twin_mNu4Sq = 92.7;
        chi2 = 'chi2Stat';
end
DataType = 'Twin';
freePar = 'E0 Norm Bkg';
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
A.exclDataStart = A.GetexclDataStart(range);
%% configure Sterile analysis object
SterileArg = {'RunAnaObj',A,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
    'nGridSteps',nGridSteps,...
    'SmartGrid','OFF',...
    'RecomputeFlag','OFF',...
    'SysEffect','all',...
    'range',range,...
    'RandMC','OFF',...
    'Twin_mNu4Sq',Twin_mNu4Sq,...
    'Twin_sin2T4',Twin_sin2T4,...
    'InterpMode','Mix'};

S = SterileAnalysis(SterileArg{:});
%%
for i=1:numel(randMC)
    
    S.RandMC= randMC(i);
    S.LoadGridFile('mNu4SqTestGrid',2,'ExtmNu4Sq','ON');
    S.Interp1Grid;
    S.FindBestFit;
    
    % get randomized spectrum
    filename = S.GridFilename('mNu4SqTestGrid',2,'ExtmNu4Sq','ON');
    d = importdata(filename);
    A.RunData.TBDIS = d.TBDIS_mc;
    if strcmp(Hypothesis,'H0')
        FitResults_Null = d.FitResults_NullOld;
    end
    
    % null hypothesis
    A.ModelObj.SetFitBias(0);
    A.ModelObj.SetFitBiasSterile(Twin_mNu4Sq,Twin_sin2T4);
    A.Fit;
    ReCalc_chi2Null = A.FitResult.chi2min;
    
    % best fit
    A.ModelObj.SetFitBias(0);
    A.ModelObj.SetFitBiasSterile(S.mNu4Sq_bf,S.sin2T4_bf);
    A.Fit;
    ReCalc_chi2Bf = A.FitResult.chi2min;
    
    if strcmp(Hypothesis,'H0')
        FitResults_NullMed = [];
        FitResults_NullOld = [];
        save(filename,'ReCalc_chi2Null','ReCalc_chi2Bf',...
            'FitResults_Null','FitResults_NullMed','FitResults_NullOld','-append');
    else
        save(filename,'ReCalc_chi2Null','ReCalc_chi2Bf','-append');
    end
end



 
 
 

