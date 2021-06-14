% ksn2 calculate chi2 grid search
%% settings that might change
chi2 = 'chi2CMShape';
DataType = 'Real';
range = 40;
freePar = 'mNu E0 Norm Bkg';
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
    'SmartGrid','OFF',...
    'RecomputeFlag','OFF',...
    'SysEffect','all',...
    'RandMC','OFF',...
    'range',range,...
    'LoadGridArg',{'ExtmNu4Sq','ON','mNu4SqTestGrid',5}};

%%
S = SterileAnalysis(SterileArg{:});
%%
if contains(freePar,'mNu')
    S.nGridSteps = 50;
    S.LoadGridFile('ExtmNu4Sq','OFF','mNu4SqTestGrid',5,'Extsin2T4','OFF');
    S.InterpMode = 'spline';
else
    S.nGridSteps = 30;
    S.LoadGridFile('ExtmNu4Sq','ON','mNu4SqTestGrid',5,'Extsin2T4','ON');
    S.InterpMode = 'Mix';
end
S.Interp1Grid;
S.ContourPlot('BestFit','ON','SavePlot','OFF');
S.CompareBestFitNull
%S.GridPlot('Contour','ON','BestFit','ON','SavePlot','png');
%%
S.FindBestFit;
fprintf('Best fit: sin2T4 = %.2f m4Sq = %.2feV^2 mNuSq = %.2feV^2 chi2min = %.2f (%.0f dof) \n',...
    S.sin2T4_bf,S.mNu4Sq_bf,S.mNuSq_bf,S.chi2_bf,S.dof);
S.FindBestFit('Mode','Imp');
fprintf('Imp. Best fit: sin2T4 = %.2f m4Sq = %.2feV^2 mNuSq = %.2feV^2 chi2min = %.2f (%.0f dof) \n',...
    S.sin2T4_bf,S.mNu4Sq_bf,S.mNuSq_bf,S.chi2_bf,S.dof);

%% 
if ~contains(freePar,'mNu')
S.LoadGridFile('ExtmNu4Sq','ON','mNu4SqTestGrid',5,'Extsin2T4','ON');
S.InterpMode = 'Mix';
S.Interp1Grid('Maxm4Sq',36^2);
 S.ContourPlot('BestFit','ON','SavePlot','ON','ExtraStr','_Extsin2T4');
 S.GridPlot('Contour','ON','BestFit','ON','SavePlot','png','ExtraStr','_Extsin2T4');
 S.FindBestFit;
fprintf('Best fit: sin2T4 = %.2f m4Sq = %.2feV^2 mNuSq = %.2feV^2 chi2min = %.2f (%.0f dof) \n',...
    S.sin2T4_bf,S.mNu4Sq_bf,S.mNuSq_bf,S.chi2_bf,S.dof);
S.FindBestFit('Mode','Imp');
fprintf('Imp. Best fit: sin2T4 = %.2f m4Sq = %.2feV^2 mNuSq = %.2feV^2 chi2min = %.2f (%.0f dof) \n',...
    S.sin2T4_bf,S.mNu4Sq_bf,S.mNuSq_bf,S.chi2_bf,S.dof);

end
%    BF = 'OFF';
% end
% S.Interp1Grid('RecomputeFlag','ON','Maxm4Sq',33^2);
% %S.GridPlot('Contour','ON','BestFit',BF,'SavePlot','png','CL',95)
% S.FindBestFit;
% S.FindBestFit('Mode','Imp');
% 
% mNu4Sq_bf = S.mNu4Sq_bf;
% sin2T4_bf = S.sin2T4_bf;
% chi2_bf = S.chi2_bf;
% %% enlarge sin2t4
% S.LoadGridFile('ExtmNu4Sq','ON','mNu4SqTestGrid',5,'Extsin2T4','ON','IgnoreKnm2FSDbinning','ON');
% 
% if strcmp(A.DataType,'Real')
%     S.InterpMode = 'spline';
%     BF = 'ON';
% else
%     S.InterpMode = 'spline';
%    BF = 'OFF';
% end
% S.Interp1Grid('RecomputeFlag','ON','Maxm4Sq',34^2);
% S.GridPlot('Contour','ON','BestFit',BF,'SavePlot','png','CL',95)
% S.FindBestFit;
% S.FindBestFit('Mode','Imp');
% mNu4Sq_bf2 = S.mNu4Sq_bf;
% sin2T4_bf2 = S.sin2T4_bf;
% chi2_bf2 = S.chi2_bf;
