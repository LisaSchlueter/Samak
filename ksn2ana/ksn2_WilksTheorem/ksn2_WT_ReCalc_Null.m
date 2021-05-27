% recalculate chi2 for MC truth
% something weird must have happened...

% grid search on randomized twins
% ksn2 calculae chi2 grid search

Hypothesis = 'H1';
switch Hypothesis
    case 'H0'
        randMC = 1:1e3;
        Twin_sin2T4 = 0;
        Twin_mNu4Sq = 0;
        chi2 = 'chi2CMShape';
    case 'H1'
        randMC = [1:139,577:757];
        Twin_sin2T4 = 0.0240;
        Twin_mNu4Sq = 92.7;
        chi2 = 'chi2CMShape';
end
DataType    = 'Twin';
freePar     = 'E0 Norm Bkg';
nGridSteps  = 25;
range       = 40;

%% get randomized MC data
savedir = [getenv('SamakPath'),'ksn2ana/ksn2_WilksTheorem/results/'];
if Twin_sin2T4==0 && Twin_mNu4Sq==0
    savefile = sprintf('%sksn2_WilksTheorem_RandomSpectra_%s_%.0feV_NullHypothesis_%.0fsamples.mat',savedir,chi2,range,5000);
else
    savefile = sprintf('%sksn2_WilksTheorem_RandomSpectra_%s_%.0feV_mNu4Sq-%.1feV2_sin2T4-%.3g_%.0fsamples.mat',...
        savedir,chi2,range,Twin_mNu4Sq,Twin_sin2T4,5000);
end

if exist(savefile,'file')
    d = importdata(savefile);
else
    fprintf('file %s not found \n',savefile)
    return
end

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
TBDIS_i = A.RunData.TBDIS;
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
    'InterpMode','lin'};

S = SterileAnalysis(SterileArg{:});
%%
for i=1:numel(randMC)
    S.RandMC = randMC(i);
    S.RandMC_TBDIS =   d.TBDIS_mc(:,randMC(i));
    filename = S.GridFilename('mNu4SqTestGrid',2,'ExtmNu4Sq','ON');
    dtmp = importdata(filename);
    
    % re-init RunAna
    A.InitModelObj_Norm_BKG('RecomputeFlag','ON');
   % A.ComputeCM;
    A.FitCMShape = dtmp.FitCMShape;
    A.RunData.TBDIS = d.TBDIS_mc(:,randMC(i));
    
    % null hypothesis with recalculated
    A.ModelObj.SetFitBias(0);
    A.ModelObj.SetFitBiasSterile(Twin_mNu4Sq,Twin_sin2T4);  
    A.Fit;
    ReCalc_chi2Null = A.FitResult.chi2min;
    
    if strcmp(Hypothesis,'H0')
        ReCalc_chi2Bf = [];
        save(filename,'ReCalc_chi2Null','ReCalc_chi2Null_i','ReCalc_chi2Bf','-append');
    else
        save(filename,'ReCalc_chi2Null','-append');
    end
end



 
 
 

