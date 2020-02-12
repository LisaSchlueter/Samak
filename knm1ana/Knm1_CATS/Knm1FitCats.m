%% settings
RunList               = 'KNM1';
exclDataStart         = 14;
RecomputeFlag         = 'ON';
SysEffects            = struct('TASR','ON','FSD','ON','RF_RX','ON','RF_EL','ON','RF_BF','ON','BkgShape','ON','TCoff_RAD','ON','TCoff_OTHER','ON');
NonPoissonScaleFactor = 1.1;
 
%% Init Model Object and covariance matrix object
Twin = MultiRunAnalysis('RunList',RunList,...
    'chi2','chi2CMShape','DataType','Twin',...
    'exclDataStart',exclDataStart,...
    'fixPar','5 6 7 8 9 10 11',...
    'RadiativeFlag','ON',...
    'NonPoissonScaleFactor',NonPoissonScaleFactor);
%Twin.ComputeCM('SysEffects',SysEffects,'BkgCM','ON');
qUmin = round(Twin.ModelObj.qU(Twin.exclDataStart));
qUmax = round(Twin.ModelObj.qU(end));
range = round(Twin.ModelObj.qU(Twin.exclDataStart)-Twin.ModelObj.Q_i);

%% Save Asimov Data and Fluctuation Statistic Only
%Twin_TBDIS_Asimov   = Twin.RunData.TBDIS;
%Twin_TBDISE_Asimov  = Twin.RunData.TBDISE;
Twin.RunData.TBDIS  = Twin_TBDIS_Asimov + randn(numel(Twin.RunData.qU),1).*Twin_TBDISE_Asimov;
Twin.RunData.TBDISE = sqrt(Twin.RunData.TBDIS);

%% Fit 
Twin.minuitOpt='migrad;minos';
Twin.exclDataStart=exclDataStart;
Twin.fixPar=' 5 6 7 8 9 10 11';
Twin.Fit('CATS','ON');
Twin.PlotFit;

%% Display Effect
Twin.FitCM_Obj.DisplayCMInfo

%% Plot Total Covariance Matrix / Correlation MAtrix
Twin.PlotFitCovCorMatrices('Mode','Frac');
Twin.PlotFitCovCorMatrices('Mode','Shape');
Twin.PlotFitCovCorMatrices('Mode','CM');
