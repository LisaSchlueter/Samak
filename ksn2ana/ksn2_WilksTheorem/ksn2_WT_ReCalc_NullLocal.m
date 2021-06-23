% recalculate chi2 for MC truth
% something weird must have happened...

% grid search on randomized twins
% ksn2 calculae chi2 grid search

Hypothesis = 'H1';
switch Hypothesis
    case 'H0'
        Twin_sin2T4 = 0;
        Twin_mNu4Sq = 0;
        chi2 = 'chi2CMShape';
    case 'H1'
           randMC = [1:129,578:748];
        Twin_sin2T4 = 0.0240;
        Twin_mNu4Sq = 92.7;
        chi2 = 'chi2CMShape';
        MergeNew = 'OFF'; % nothing new
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
S = SterileAnalysis(SterileArg{:});

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
%% get grid
S.RandMC = 3;
S.RandMC_TBDIS = d.TBDIS_mc(:,3);
S.LoadGridFile('mNu4SqTestGrid',2,'ExtmNu4Sq','ON');
S.FindBestFit;
%%

% S.RandMC= 3;
% filename = S.GridFilename('mNu4SqTestGrid',2,'ExtmNu4Sq','ON');
% %d = importdata(filename);

%A.ComputeCM
% %
% mNu4SqTest = d.mnu4Sq(4,25);
% sin2T4Test =  d.sin2T4(4,25);
% chi2 = d.chi2';
% chi2Test = chi2(4,25);

%%
%A.RunData.TBDIS = d.TBDIS_mc;
A.RunData.TBDIS = [zeros(10,1)',d.TBDIS_i]';
A.ComputeCM
A.RunData.TBDIS = d.TBDIS_mc(:,3);
A.ModelObj.SetFitBiasSterile(Twin_mNu4Sq,Twin_sin2T4);%mNu4SqTest,sin2T4Test);
A.Fit;
chi2_ReCalc = A.FitResult.chi2min;
fprintf('chi2null from grid search  : %.2f \n',S.chi2_Null) 
fprintf('chi2null from recalculating: %.2f \n',chi2_ReCalc) 
%%
A.RunData.TBDIS = TBDIS_i;
A.ComputeCM;
A.ModelObj.SetFitBiasSterile(1,0.5);%mNu4SqTest,sin2T4Test);
A.RunData.TBDIS = d.TBDIS_mc;
A.Fit;
chi2_ReCalc2 = A.FitResult.chi2min;


% grid that have best fit at null
a = [  3
     8
    17
    23
    26
    28
    48
    66
    95
   119
   127
   137
   144
   158
   170
   207
   211
   216
   227
   237
   244
   258
   270
   303
   317
   321
   323
   351
   361
   370
   375
   376
   383
   388
   395
   403
   410
   432
   433
   434
   443
   444
   449
   453
   470
   502
   507
   514
   522
   529
   541
   543
   558
   560
   566
   573
   574
   582
   601
   605
   606
   609
   612
   628
   637
   644
   658
   670
   703
   705
   713
   727
   734
   738
   739
   748
   758
   767
   769
   774
   784
   788
   792
   794
   804
   810
   812
   813
   830
   841
   849
   856
   864
   902
   914
   928
   935
   947
   962
   975];
