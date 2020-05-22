%% cross check best fit again
chi2Str = 'chi2CMShape';
range = 65;
sin_bf =  0.0251;
m_bf= 78.1921;

%% martins's best fit
mM_bf   = 201.984663;
sinM_bf = 0.008273;
%% 1. RunAnalysis object
RunAnaArg = {'RunList','KNM1',...
    'fixPar','E0 Norm Bkg',...
    'DataType','Real',...
    'FSDFlag','SibilleFull',...
    'ELossFlag','KatrinT2',...
    'AnaFlag','StackPixel',...
    'chi2','chi2Stat',...
    'ROIFlag','Default',...
    'SynchrotronFlag','ON',...
    'AngularTFFlag','OFF',...
    'ISCSFlag','Edep',...
    'TwinBias_Q',18573.73,...
    'SysBudget',24,...
    'pullFlag',99,...
    'NonPoissonScaleFactor',1};

R = MultiRunAnalysis(RunAnaArg{:});
R.exclDataStart = R.GetexclDataStart(range);
R.chi2 = chi2Str;

if strcmp(chi2Str,'chi2CMShape')
    R.NonPoissonScaleFactor=1.064;
    R.ComputeCM;
else
    R.NonPoissonScaleFactor=1;
end

%% samak best fit
R.ModelObj.SetFitBiasSterile(m_bf,sin_bf);
R.Fit;
FitResultS = R.FitResult;
%% martin best fit
R.ModelObj.SetFitBiasSterile(mM_bf,sinM_bf);
R.Fit;
FitResultM = R.FitResult;

%% compare chi2
fprintf('----------------------\n');
fprintf('chi2min Samak   %.2f \n',FitResultS.chi2min);
fprintf('chi2min Fitrium %.2f \n',FitResultM.chi2min);
fprintf('----------------------\n');
