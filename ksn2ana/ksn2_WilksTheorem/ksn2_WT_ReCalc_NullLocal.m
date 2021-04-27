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
%% configure Sterile analysis object
SterileArg = {'RunAnaObj',A,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
    'nGridSteps',nGridSteps,...
    'SmartGrid','OFF',...
    'RecomputeFlag','OFF',...
    'SysEffect','all',...
    'range',range,...
    'RandMC','OFF',...
    'Twin_mNu4Sq',Twin_mNu4Sq,...
    'Twin_sin2T4',Twin_sin2T4};
A = MultiRunAnalysis(RunAnaArg{:});
A.exclDataStart = A.GetexclDataStart(range);
A.InitModelObj_Norm_BKG('RecomputeFlag','ON');
TBDIS_i = A.RunData.TBDIS;
%%

% index of grid (h0) that have delta chi^2 > 5.99, but are outside of contour (sensitivity)
myRandMC = [  1
    17
    48
    63
    82
    132
    165
    174
    176
    198
    200
    232
    265
    274
    276
    298
    300
    307
    375
    381
    422
    435
    444
    446
    474
    476
    498
    500
    573
    590
    618
    632
    665
    674
    676
    698
    700
    719
    759
    760
    802
    864
    871
    874
    896
    898
    907
    940
    966
    970
    980
    982
    988
    1000 ];

chi2NullFile = zeros(numel(myRandMC),1);
chi2BfFile = zeros(numel(myRandMC),1);
chi2NullReCalc = zeros(numel(myRandMC),1);

for j=1:1:numel(myRandMC)
    progressbar(j/numel(myRandMC));
   
    % load grid and find best fit & null hypothesis
    S = SterileAnalysis(SterileArg{:});
    S.RandMC= myRandMC(j);%randMC(i);
    S.InterpMode = 'lin';
    S.LoadGridFile('mNu4SqTestGrid',2,'ExtmNu4Sq','ON');
    S.Interp1Grid('Minm4Sq',1);
    S.FindBestFit;
    chi2NullFile(j)   = S.chi2_Null;
    chi2BfFile(j)     = S.chi2_bf;
    
 % load file
    filename = S.GridFilename('mNu4SqTestGrid',2,'ExtmNu4Sq','ON');
    d = importdata(filename);
      
    %% check null hypothesis with init from TBDIS_i
    A.RunData.TBDIS = TBDIS_i;
    A.ModelObj.SetFitBiasSterile(0,0);
    A.InitModelObj_Norm_BKG('RecomputeFlag','OFF');
    
    A.RunData.TBDIS = d.TBDIS_mc;
    A.ModelObj.SetFitBiasSterile(Twin_mNu4Sq,Twin_sin2T4);
    
    A.Fit;
    chi2NullReCalc(j) = A.FitResult.chi2min;
    
end

%%
chi2deltaFile = chi2NullFile - chi2BfFile;
SigFracFile = 1e2*sum(chi2deltaFile>=5.99)./numel(chi2deltaFile);

chi2deltaNew = chi2NullReCalc - chi2BfFile;
SigFracNew =1e2*sum(chi2deltaNew>=5.99)./numel(chi2deltaFile);

%% cross check chi2min from best fit
% A.ModelObj.SetFitBiasSterile(S.mNu4Sq_bf,S.sin2T4_bf);
% A.Fit;
% chi2BfFile   = S.chi2_bf;
% chi2BfReCalc = A.FitResult.chi2min;
% 
% %%
% FitResults_NullOld = d.FitResults_Null;
% FitResults_Null = A.FitResult;
% 
% 
