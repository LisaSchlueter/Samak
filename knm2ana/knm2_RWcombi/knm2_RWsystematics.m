% Test new rear wall potential (RW) covariance matrix

%% set up model object
RunList = 'KNM2_Prompt';
FSDFlag = 'BlindingKNM2';
ELossFlag = 'KatrinT2';
AnaFlag = 'StackPixel'; % uniform FPD
chi2 = 'chi2Stat';
DataType = 'Twin';
freePar = 'mNu E0 Bkg Norm';
range = 40;
RunArg =     {'FSDFlag',FSDFlag,...
    'ELossFlag',ELossFlag,...
    'AnaFlag',AnaFlag,...
    'chi2',chi2,...
    'RunList',RunList,...
    'fixPar',freePar};  % twins, all runs have same endpoint

T  = MultiRunAnalysis(RunArg{:},'DataType','Twin','TwinBias_Q','Fit'); % Twins with gaussian variation of endpoint
T.exclDataStart =  T.GetexclDataStart(range);
RunArg = {RunArg{:},'exclDataStart',T.exclDataStart};

%% covariance matrix object
CMarg = {'StudyObject',T.ModelObj,...
         'nTrials',100,...
         'SysEffect',struct('RW','ON'),...
         'RecomputeFlag','ON',...
         'SanityPlots','OFF',...
         'RW_SigmaErr',0.1};

CM = CovarianceMatrix(CMarg{:});
CM.ComputeCM_RW('Sigma',0.2,'Dist','Gauss');
%%
CM.PlotCM('qUWindowIndexMin',90,'qUWindowIndexMax',40,'Convergence','ON')